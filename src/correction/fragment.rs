//! Fragment merging and mistranslation correction.
//!
//! Compares nearby genes at the nucleotide level and merges
//! clusters with identical DNA sequences (>=95% coverage, >=99% identity).
//! Also collapses gene families sharing common neighbors at configurable
//! thresholds, iterating from high to low identity (matching Panaroo's
//! progressive collapsing behavior).

use std::collections::{HashMap, HashSet};

use crate::error::Result;
use crate::graph::ConcurrentGraph;
use crate::graph::ClusterId;
use crate::correction::simd::align_sequences;

/// Cached pairwise distances between cluster centroids.
///
/// Stores computed identity scores so they can be reused across
/// correction passes (e.g., after missing gene recovery), matching
/// Panaroo's approach of reusing the distance matrix.
#[derive(Debug, Clone, Default)]
pub struct DistanceCache {
    /// Maps (cluster_a, cluster_b) → identity score
    /// Keys are normalized so a < b lexicographically
    distances: HashMap<(String, String), f64>,
}

impl DistanceCache {
    /// Create a new empty distance cache.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a distance between two clusters (order-independent key).
    pub fn insert(&mut self, a: &str, b: &str, identity: f64) {
        let key = Self::normalize_key(a, b);
        self.distances.insert(key, identity);
    }

    /// Look up a cached distance between two clusters.
    /// Avoids allocation by scanning entries with &str comparison.
    pub fn get(&self, a: &str, b: &str) -> Option<f64> {
        let (lo, hi) = Self::normalize_key_parts(a, b);
        self.distances.iter()
            .find(|((k1, k2), _)| k1.as_str() == lo && k2.as_str() == hi)
            .map(|(_, v)| *v)
    }

    /// Get the number of cached distances.
    pub fn len(&self) -> usize {
        self.distances.len()
    }

    /// Check if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.distances.is_empty()
    }

    /// Normalize key so the smaller string is first (order-independent lookup).
    fn normalize_key(a: &str, b: &str) -> (String, String) {
        if a < b {
            (a.to_string(), b.to_string())
        } else {
            (b.to_string(), a.to_string())
        }
    }

    /// Returns (lo, hi) without allocating, for lookup.
    fn normalize_key_parts<'a>(a: &'a str, b: &'a str) -> (&'a str, &'a str) {
        if a < b { (a, b) } else { (b, a) }
    }
}

/// Merges fragmented genes and corrects mistranslations.
pub struct FragmentMerger {
    /// Minimum coverage for mistranslation correction (default: 0.95)
    coverage_threshold: f32,
    /// Minimum identity for mistranslation correction (default: 0.99)
    identity_threshold: f32,
    /// Thresholds for iterative gene family collapsing, high to low (default: [0.99, 0.95, 0.9, 0.8, 0.7])
    collapse_thresholds: Vec<f32>,
    /// BFS depth for neighbor search (default: 3, matching Panaroo)
    bfs_depth: usize,
}

impl FragmentMerger {
    /// Create a new fragment merger with default thresholds.
    pub fn new() -> Self {
        Self {
            coverage_threshold: 0.95,
            identity_threshold: 0.99,
            collapse_thresholds: vec![0.99, 0.95, 0.9, 0.8, 0.7],
            bfs_depth: 3,
        }
    }

    /// Set custom collapsing thresholds (high to low).
    pub fn with_collapse_thresholds(mut self, thresholds: Vec<f32>) -> Self {
        self.collapse_thresholds = thresholds;
        self
    }

    /// Backward-compatible: set a single collapsing threshold.
    pub fn with_collapse_threshold(mut self, threshold: f32) -> Self {
        self.collapse_thresholds = vec![threshold];
        self
    }

    /// Set BFS depth for neighbor search.
    pub fn with_bfs_depth(mut self, depth: usize) -> Self {
        self.bfs_depth = depth;
        self
    }

    /// Get the configured collapsing thresholds.
    pub fn collapse_thresholds(&self) -> &[f32] {
        &self.collapse_thresholds
    }

    /// Run mistranslation correction on the graph with BFS neighbor search.
    ///
    /// Searches at BFS depths [1..bfs_depth] for similar clusters,
    /// matching Panaroo's approach of checking neighbors at multiple depths.
    /// Uses Levenshtein alignment instead of positional comparison to
    /// correctly handle insertions and deletions.
    pub fn correct_mistranslations(
        &self,
        graph: &ConcurrentGraph,
        sequences: &HashMap<String, Vec<u8>>,
    ) -> Result<MergeStats> {
        let mut merged_count = 0;

        // Collect all cluster IDs
        let cluster_ids: Vec<ClusterId> = graph.nodes.iter()
            .map(|e| e.key().clone())
            .collect();

        for cluster_id in &cluster_ids {
            // Search neighbors at configured BFS depth (matching Panaroo's [1,2,3])
            let neighbors = self.collect_neighbors_at_depth(graph, cluster_id, self.bfs_depth);

            for neighbor_id in &neighbors {
                // Skip if either node was already merged
                if graph.nodes.get(cluster_id).is_none() || graph.nodes.get(neighbor_id).is_none() {
                    continue;
                }

                let seq_a = sequences.get(cluster_id.as_str());
                let seq_b = sequences.get(neighbor_id.as_str());

                if let (Some(a), Some(b)) = (seq_a, seq_b) {
                    let (coverage, identity) = Self::compare_sequences(a, b);

                    if coverage >= self.coverage_threshold && identity >= self.identity_threshold {
                        let support_a = graph.nodes.get(cluster_id).map(|n| n.support).unwrap_or(0);
                        let support_b = graph.nodes.get(neighbor_id).map(|n| n.support).unwrap_or(0);

                        if support_a >= support_b {
                            graph.merge_nodes(cluster_id, neighbor_id);
                        } else {
                            graph.merge_nodes(neighbor_id, cluster_id);
                        }
                        merged_count += 1;
                    }
                }
            }
        }

        tracing::info!("Mistranslation correction: merged {} cluster pairs", merged_count);

        Ok(MergeStats {
            pairs_merged: merged_count,
            families_collapsed: 0,
        })
    }

    /// Run gene family collapsing on the graph.
    ///
    /// Uses pre-clustering by length and prefix for efficiency (instead of O(n²)
    /// all-pairs comparison). Only compares clusters within the same length/prefix
    /// bucket that share a common neighbor.
    pub fn collapse_gene_families(
        &self,
        graph: &ConcurrentGraph,
        sequences: &HashMap<String, Vec<u8>>,
    ) -> Result<usize> {
        self.collapse_gene_families_with_cache(graph, sequences, None)
    }

    /// Run gene family collapsing with an optional distance cache.
    ///
    /// When a cache is provided, computed distances are stored for reuse
    /// in subsequent passes (matching Panaroo's distance matrix reuse).
    /// Uses the first threshold from `collapse_thresholds` for backward
    /// compatibility with single-threshold callers.
    pub fn collapse_gene_families_with_cache(
        &self,
        graph: &ConcurrentGraph,
        sequences: &HashMap<String, Vec<u8>>,
        cache: Option<&mut DistanceCache>,
    ) -> Result<usize> {
        let threshold = self.collapse_thresholds.first().copied().unwrap_or(0.7);
        self.collapse_gene_families_with_threshold(graph, sequences, threshold, cache)
    }

    /// Collapse gene families at a specific threshold.
    ///
    /// This is the core collapsing routine. Callers (e.g., the pipeline)
    /// iterate over multiple thresholds from high to low, reusing the
    /// same `DistanceCache` across iterations.
    pub fn collapse_gene_families_with_threshold(
        &self,
        graph: &ConcurrentGraph,
        sequences: &HashMap<String, Vec<u8>>,
        threshold: f32,
        mut cache: Option<&mut DistanceCache>,
    ) -> Result<usize> {
        let mut collapsed = 0;

        // Collect cluster IDs
        let cluster_ids: Vec<ClusterId> = graph.nodes.iter()
            .map(|e| e.key().clone())
            .collect();

        // Pre-cluster for efficiency (instead of O(n²) all-pairs)
        let groups = Self::pre_cluster_by_similarity(&cluster_ids, sequences, 100);

        for group in &groups {
            for i in 0..group.len() {
                for j in (i + 1)..group.len() {
                    let id_a = &cluster_ids[group[i]];
                    let id_b = &cluster_ids[group[j]];

                    // Skip if either node was already merged
                    if graph.nodes.get(id_a).is_none() || graph.nodes.get(id_b).is_none() {
                        continue;
                    }

                    // Check if they share a neighbor (skip if not connected)
                    if !self.share_neighbor(graph, id_a, id_b) {
                        continue;
                    }

                    // Check cache first
                    let identity = match &mut cache {
                        Some(c) => {
                            if let Some(cached_id) = c.get(id_a.as_str(), id_b.as_str()) {
                                cached_id as f32
                            } else {
                                let (cov, id) = Self::compare_sequences_from_seqs(
                                    id_a, id_b, sequences
                                );
                                if cov > 0.0 {
                                    c.insert(id_a.as_str(), id_b.as_str(), id as f64);
                                }
                                id
                            }
                        }
                        None => {
                            let (_, id) = Self::compare_sequences_from_seqs(
                                id_a, id_b, sequences
                            );
                            id
                        }
                    };

                    if identity >= threshold {
                        // Merge the node with lower support into the one with higher support
                        let support_a = graph.nodes.get(id_a).map(|n| n.support).unwrap_or(0);
                        let support_b = graph.nodes.get(id_b).map(|n| n.support).unwrap_or(0);

                        if support_a >= support_b {
                            graph.merge_nodes(id_a, id_b);
                        } else {
                            graph.merge_nodes(id_b, id_a);
                        }
                        collapsed += 1;
                    }
                }
            }
        }

        tracing::info!("Gene family collapsing (threshold={:.2}): collapsed {} clusters", threshold, collapsed);
        Ok(collapsed)
    }

    /// Compare two sequences and return (coverage, identity).
    ///
    /// Uses Levenshtein alignment (edit distance) instead of positional comparison,
    /// correctly handling insertions and deletions.
    fn compare_sequences(a: &[u8], b: &[u8]) -> (f32, f32) {
        if a.is_empty() || b.is_empty() {
            return (0.0, 0.0);
        }

        let min_len = a.len().min(b.len());
        let max_len = a.len().max(b.len());

        let coverage = min_len as f32 / max_len as f32;

        // Use Levenshtein alignment instead of positional comparison
        let identity = align_sequences(a, b) as f32;

        (coverage, identity)
    }

    /// Compare sequences by cluster IDs using Levenshtein alignment.
    fn compare_sequences_from_seqs(
        id_a: &ClusterId,
        id_b: &ClusterId,
        sequences: &HashMap<String, Vec<u8>>,
    ) -> (f32, f32) {
        let seq_a = sequences.get(id_a.as_str());
        let seq_b = sequences.get(id_b.as_str());

        match (seq_a, seq_b) {
            (Some(a), Some(b)) => Self::compare_sequences(a, b),
            _ => (0.0, 0.0),
        }
    }

    /// Check if two clusters share a common neighbor in the graph.
    /// Uses the adjacency index for O(degree) lookups instead of O(E) scans.
    fn share_neighbor(
        &self,
        graph: &ConcurrentGraph,
        id_a: &ClusterId,
        id_b: &ClusterId,
    ) -> bool {
        let neighbors_a: HashSet<ClusterId> = graph.neighbors(id_a).into_iter().collect();
        let neighbors_b = graph.neighbors(id_b);
        neighbors_b.iter().any(|n| neighbors_a.contains(n))
    }

    /// Pre-cluster gene families by length and prefix for efficient comparison.
    ///
    /// Instead of O(n²) all-pairs comparison, groups clusters into buckets
    /// by (length_bucket, 4-byte prefix) so only similar sequences are compared.
    fn pre_cluster_by_similarity(
        cluster_ids: &[ClusterId],
        sequences: &HashMap<String, Vec<u8>>,
        bucket_size: usize,
    ) -> Vec<Vec<usize>> {
        let mut buckets: HashMap<(usize, [u8; 4]), Vec<usize>> = HashMap::new();

        for (idx, id) in cluster_ids.iter().enumerate() {
            if let Some(seq) = sequences.get(id.as_str()) {
                let len_bucket = seq.len() / bucket_size;
                let prefix = if seq.len() >= 4 {
                    [seq[0], seq[1], seq[2], seq[3]]
                } else if !seq.is_empty() {
                    let mut p = [0u8; 4];
                    for (i, &b) in seq.iter().take(4).enumerate() {
                        p[i] = b;
                    }
                    p
                } else {
                    [0u8; 4]
                };
                buckets.entry((len_bucket, prefix)).or_default().push(idx);
            }
        }

        buckets.into_values().collect()
    }

    /// Collect neighbors within a given BFS depth.
    ///
    /// This implements Panaroo's approach of searching at depths [1,2,3]
    /// for similar clusters, instead of only checking edge-adjacent pairs.
    /// Uses the adjacency index for O(degree) lookups instead of O(E) scans.
    fn collect_neighbors_at_depth(
        &self,
        graph: &ConcurrentGraph,
        start: &ClusterId,
        max_depth: usize,
    ) -> Vec<ClusterId> {
        let mut visited = std::collections::HashSet::new();
        visited.insert(start.clone());

        let mut current_level = vec![start.clone()];
        let mut neighbors = Vec::new();

        for _depth in 0..max_depth {
            let mut next_level = Vec::new();
            for node_id in &current_level {
                // Use adjacency index for O(degree) lookup instead of O(E) scan
                for neighbor in graph.neighbors(node_id) {
                    if !visited.contains(&neighbor) {
                        visited.insert(neighbor.clone());
                        next_level.push(neighbor.clone());
                        neighbors.push(neighbor);
                    }
                }
            }
            current_level = next_level;
        }

        neighbors
    }
}

impl Default for FragmentMerger {
    fn default() -> Self {
        Self::new()
    }
}

/// Statistics from fragment merging.
#[derive(Debug, Clone)]
pub struct MergeStats {
    /// Number of cluster pairs merged (mistranslation correction)
    pub pairs_merged: usize,
    /// Number of gene families collapsed
    pub families_collapsed: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compare_sequences() {
        let (cov, id) = FragmentMerger::compare_sequences(b"ATCGATCG", b"ATCGATCG");
        assert_eq!(cov, 1.0);
        // Levenshtein alignment gives identity 1.0 for identical sequences
        assert!((id - 1.0).abs() < 0.01, "Expected ~1.0, got {}", id);

        let (cov, id) = FragmentMerger::compare_sequences(b"ATCGATCG", b"ATCGATCC");
        assert_eq!(cov, 1.0);
        // Levenshtein: 1 substitution in 8 chars = 7/8 = 0.875
        assert!((id - 0.875).abs() < 0.01, "Expected ~0.875, got {}", id);
    }

    #[test]
    fn test_compare_empty() {
        let (cov, id) = FragmentMerger::compare_sequences(b"", b"ATCG");
        assert_eq!(cov, 0.0);
        assert_eq!(id, 0.0);
    }

    #[test]
    fn test_distance_cache() {
        let mut cache = DistanceCache::new();
        assert!(cache.is_empty());

        cache.insert("c1", "c2", 0.85);
        assert_eq!(cache.len(), 1);

        // Order-independent lookup
        assert_eq!(cache.get("c1", "c2"), Some(0.85));
        assert_eq!(cache.get("c2", "c1"), Some(0.85));

        // Non-existent entry
        assert_eq!(cache.get("c1", "c3"), None);
    }

    #[test]
    fn test_distance_cache_order_independence() {
        let mut cache = DistanceCache::new();
        cache.insert("b", "a", 0.99);

        // Both orderings should return the same value
        assert_eq!(cache.get("a", "b"), Some(0.99));
        assert_eq!(cache.get("b", "a"), Some(0.99));
    }

    #[test]
    fn test_pre_cluster_by_similarity() {
        let mut sequences = HashMap::new();
        sequences.insert("c1".to_string(), b"ATCGATCGATCGATCG".to_vec()); // length 16
        sequences.insert("c2".to_string(), b"ATCGATCGATCGATCC".to_vec()); // length 16, same prefix
        sequences.insert("c3".to_string(), b"GGGGGGGGGGGGGGGG".to_vec()); // length 16, different prefix
        sequences.insert("c4".to_string(), b"A".to_vec()); // length 1

        let cluster_ids: Vec<ClusterId> = vec![
            ClusterId::new("c1"),
            ClusterId::new("c2"),
            ClusterId::new("c3"),
            ClusterId::new("c4"),
        ];

        let groups = FragmentMerger::pre_cluster_by_similarity(&cluster_ids, &sequences, 100);

        // c1 and c2 should be in the same bucket (same length bucket, same prefix)
        // c3 should be in a different bucket (different prefix)
        // c4 should be in its own bucket (different length)
        assert!(groups.len() >= 2, "Expected at least 2 groups, got {}", groups.len());
    }

    #[test]
    fn test_collect_neighbors_at_depth() {
        use crate::graph::{Node, GeneCluster, GenomeId};

        let graph = ConcurrentGraph::with_capacity(10);

        // Build a chain: a -> b -> c -> d
        let nodes = vec![
            Node::from_cluster(&{ let mut c = GeneCluster::new("a"); c.support = 5; c }),
            Node::from_cluster(&{ let mut c = GeneCluster::new("b"); c.support = 5; c }),
            Node::from_cluster(&{ let mut c = GeneCluster::new("c"); c.support = 5; c }),
            Node::from_cluster(&{ let mut c = GeneCluster::new("d"); c.support = 5; c }),
        ];
        for node in nodes {
            graph.add_node(node);
        }

        graph.add_edge_genome(ClusterId::new("a"), ClusterId::new("b"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("b"), ClusterId::new("c"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("c"), ClusterId::new("d"), GenomeId::new("g1"));

        let merger = FragmentMerger::new();

        // Depth 1: a's neighbors are {b}
        let depth1 = merger.collect_neighbors_at_depth(&graph, &ClusterId::new("a"), 1);
        assert_eq!(depth1.len(), 1);

        // Depth 2: a's neighbors are {b, c}
        let depth2 = merger.collect_neighbors_at_depth(&graph, &ClusterId::new("a"), 2);
        assert_eq!(depth2.len(), 2);

        // Depth 3: a's neighbors are {b, c, d}
        let depth3 = merger.collect_neighbors_at_depth(&graph, &ClusterId::new("a"), 3);
        assert_eq!(depth3.len(), 3);
    }

    #[test]
    fn test_iterative_collapse_merges_progressively() {
        use crate::graph::{Node, GeneCluster, GenomeId};

        let graph = ConcurrentGraph::with_capacity(10);

        // Create three clusters: c1 and c2 have identical sequences and
        // share a common neighbor (c3). The collapse method requires clusters
        // to share a common neighbor before merging.
        let seq = b"ATCGATCGATCGATCGATCGATCGATCGATCG".to_vec(); // 32 bases
        let seq3 = b"GGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGG".to_vec(); // different sequence for c3

        let nodes = vec![
            Node::from_cluster(&{ let mut c = GeneCluster::new("c1"); c.support = 10; c }),
            Node::from_cluster(&{ let mut c = GeneCluster::new("c2"); c.support = 5; c }),
            Node::from_cluster(&{ let mut c = GeneCluster::new("c3"); c.support = 8; c }),
        ];
        for node in nodes {
            graph.add_node(node);
        }

        // c3 is a common neighbor of both c1 and c2
        graph.add_edge_genome(ClusterId::new("c1"), ClusterId::new("c3"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("c2"), ClusterId::new("c3"), GenomeId::new("g1"));

        let mut sequences = HashMap::new();
        sequences.insert("c1".to_string(), seq.clone());
        sequences.insert("c2".to_string(), seq);
        sequences.insert("c3".to_string(), seq3);

        let merger = FragmentMerger::new()
            .with_collapse_thresholds(vec![0.99, 0.95, 0.9, 0.8, 0.7]);

        // At threshold 0.99, identical sequences (c1+c2) should merge
        let mut cache = DistanceCache::new();
        let merged_099 = merger.collapse_gene_families_with_threshold(
            &graph, &sequences, 0.99, Some(&mut cache)
        ).unwrap();
        assert!(merged_099 >= 1, "At least 1 merge at 0.99 (identical sequences)");

        // The cache should have entries after collapsing
        assert!(!cache.is_empty(), "DistanceCache should have entries after collapsing");
    }

    #[test]
    fn test_with_collapse_threshold_backward_compat() {
        let merger = FragmentMerger::new().with_collapse_threshold(0.5);
        assert_eq!(merger.collapse_thresholds(), &[0.5]);
    }

    #[test]
    fn test_with_collapse_thresholds_multi() {
        let merger = FragmentMerger::new()
            .with_collapse_thresholds(vec![0.99, 0.95, 0.9, 0.8, 0.7]);
        assert_eq!(merger.collapse_thresholds(), &[0.99, 0.95, 0.9, 0.8, 0.7]);
    }

    #[test]
    fn test_default_collapse_thresholds() {
        let merger = FragmentMerger::new();
        assert_eq!(merger.collapse_thresholds(), &[0.99, 0.95, 0.9, 0.8, 0.7]);
    }
}