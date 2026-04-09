//! Fragment merging and mistranslation correction.
//!
//! Compares nearby genes at the nucleotide level and merges
//! clusters with identical DNA sequences (>=95% coverage, >=99% identity).
//! Also collapses gene families sharing common neighbors at 70% threshold.

use crate::error::Result;
use crate::graph::ConcurrentGraph;
use crate::graph::ClusterId;

/// Merges fragmented genes and corrects mistranslations.
pub struct FragmentMerger {
    /// Minimum coverage for mistranslation correction (default: 0.95)
    coverage_threshold: f32,
    /// Minimum identity for mistranslation correction (default: 0.99)
    identity_threshold: f32,
    /// Threshold for gene family collapsing (default: 0.70)
    collapse_threshold: f32,
}

impl FragmentMerger {
    /// Create a new fragment merger with default thresholds.
    pub fn new() -> Self {
        Self {
            coverage_threshold: 0.95,
            identity_threshold: 0.99,
            collapse_threshold: 0.70,
        }
    }

    /// Set the collapse threshold for gene family merging.
    pub fn with_collapse_threshold(mut self, threshold: f32) -> Self {
        self.collapse_threshold = threshold;
        self
    }

    /// Run mistranslation correction on the graph.
    ///
    /// Finds pairs of clusters with highly similar sequences
    /// and merges them into a single cluster.
    pub fn correct_mistranslations(
        &self,
        graph: &ConcurrentGraph,
        sequences: &std::collections::HashMap<String, Vec<u8>>,
    ) -> Result<MergeStats> {
        let mut merged_count = 0;

        // Collect node pairs that share edges
        let mut node_pairs = Vec::new();
        for entry in graph.edges.iter() {
            let (from, to) = entry.key();
            node_pairs.push((from.clone(), to.clone()));
        }

        // Check each pair for sequence similarity
        for (id_a, id_b) in &node_pairs {
            let seq_a = sequences.get(id_a.as_str());
            let seq_b = sequences.get(id_b.as_str());

            if let (Some(a), Some(b)) = (seq_a, seq_b) {
                let (coverage, identity) = Self::compare_sequences(a, b);

                if coverage >= self.coverage_threshold && identity >= self.identity_threshold {
                    // Merge b into a (keep higher-support node)
                    let support_a = graph.nodes.get(id_a).map(|n| n.support).unwrap_or(0);
                    let support_b = graph.nodes.get(id_b).map(|n| n.support).unwrap_or(0);

                    if support_a >= support_b {
                        graph.merge_nodes(id_a, id_b);
                    } else {
                        graph.merge_nodes(id_b, id_a);
                    }
                    merged_count += 1;
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
    /// Clusters sharing common neighbors are compared at a relaxed
    /// identity threshold and merged if similar enough.
    pub fn collapse_gene_families(
        &self,
        graph: &ConcurrentGraph,
        sequences: &std::collections::HashMap<String, Vec<u8>>,
    ) -> Result<usize> {
        let mut collapsed = 0;

        // Find clusters sharing common neighbors
        let mut cluster_ids = Vec::new();
        for entry in graph.nodes.iter() {
            cluster_ids.push(entry.key().clone());
        }

        for i in 0..cluster_ids.len() {
            for j in (i + 1)..cluster_ids.len() {
                let id_a = &cluster_ids[i];
                let id_b = &cluster_ids[j];

                // Check if they share a neighbor
                if !self.share_neighbor(graph, id_a, id_b) {
                    continue;
                }

                // Compare sequences at relaxed threshold
                if let (Some(a), Some(b)) = (
                    sequences.get(id_a.as_str()),
                    sequences.get(id_b.as_str()),
                ) {
                    let (_, identity) = Self::compare_sequences(a, b);
                    if identity >= self.collapse_threshold {
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

        tracing::info!("Gene family collapsing: collapsed {} clusters", collapsed);
        Ok(collapsed)
    }

    /// Compare two sequences and return (coverage, identity).
    fn compare_sequences(a: &[u8], b: &[u8]) -> (f32, f32) {
        if a.is_empty() || b.is_empty() {
            return (0.0, 0.0);
        }

        let min_len = a.len().min(b.len());
        let max_len = a.len().max(b.len());

        let coverage = min_len as f32 / max_len as f32;

        let matches = a.iter()
            .zip(b.iter())
            .filter(|(x, y)| x == y)
            .count();

        let identity = matches as f32 / min_len as f32;

        (coverage, identity)
    }

    /// Check if two clusters share a common neighbor in the graph.
    fn share_neighbor(
        &self,
        graph: &ConcurrentGraph,
        id_a: &ClusterId,
        id_b: &ClusterId,
    ) -> bool {
        let mut neighbors_a = Vec::new();
        for entry in graph.edges.iter() {
            let (from, to) = entry.key();
            if from == id_a {
                neighbors_a.push(to.clone());
            } else if to == id_a {
                neighbors_a.push(from.clone());
            }
        }

        // Check if any neighbor of id_a is also a neighbor of id_b
        graph.edges
            .iter()
            .any(|entry| {
                let (from, to) = entry.key();
                let neighbor = if from == id_b {
                    Some(to)
                } else if to == id_b {
                    Some(from)
                } else {
                    None
                };
                neighbor.map(|n| neighbors_a.contains(n)).unwrap_or(false)
            })
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
        assert_eq!(id, 1.0);

        let (cov, id) = FragmentMerger::compare_sequences(b"ATCGATCG", b"ATCGATCC");
        assert_eq!(cov, 1.0);
        assert_eq!(id, 0.875); // 7/8
    }

    #[test]
    fn test_compare_empty() {
        let (cov, id) = FragmentMerger::compare_sequences(b"", b"ATCG");
        assert_eq!(cov, 0.0);
        assert_eq!(id, 0.0);
    }
}
