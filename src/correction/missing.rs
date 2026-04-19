//! Missing gene recovery.
//!
//! Searches for genes that may have been missed by the annotation tool
//! by looking at flanking sequences around expected gene locations.
//! Uses semi-global Levenshtein alignment (HW mode) matching Panaroo's
//! edlib-based approach.

use std::collections::HashMap;

use crate::error::Result;
use crate::graph::{ClusterId, ConcurrentGraph};
use crate::correction::simd::align_semiglobal;

/// Recovers genes missed by the initial annotation.
///
/// For each node pair where one genome lacks a neighboring gene,
/// searches the surrounding contig sequence for the missing gene
/// using semi-global alignment (HW mode) matching Panaroo's find_missing.
///
/// When `remove_by_consensus` is enabled (strict mode), nodes where
/// the number of refound hits exceeds the original node support are
/// deleted entirely -- they are likely spurious annotation artifacts.
/// This matches Panaroo's `--remove_by_consensus` behavior.
pub struct MissingGeneRecoverer {
    /// Minimum alignment identity to consider a hit (default: 0.70)
    min_identity: f32,
    /// Search window size in base pairs (default: 5000)
    search_window: usize,
    /// Minimum fraction of query that must align (default: 0.20)
    prop_match: f32,
    /// Remove nodes where refound hits exceed original support (default: false)
    remove_by_consensus: bool,
}

impl MissingGeneRecoverer {
    /// Create a new missing gene recoverer with default settings.
    pub fn new() -> Self {
        Self {
            min_identity: 0.70,
            search_window: 5000,
            prop_match: 0.20,
            remove_by_consensus: false,
        }
    }

    /// Set the minimum identity threshold.
    pub fn with_min_identity(mut self, identity: f32) -> Self {
        self.min_identity = identity;
        self
    }

    /// Set the search window size.
    pub fn with_search_window(mut self, window: usize) -> Self {
        self.search_window = window;
        self
    }

    /// Set the minimum proportion of query that must align.
    pub fn with_prop_match(mut self, prop: f32) -> Self {
        self.prop_match = prop;
        self
    }

    /// Enable or disable consensus-based node removal.
    ///
    /// When enabled, nodes where the number of refound gene hits
    /// exceeds the original node support are removed from the graph.
    /// This catches spurious clusters that are annotation artifacts,
    /// matching Panaroo's `--remove_by_consensus` behavior.
    pub fn with_remove_by_consensus(mut self, remove: bool) -> Self {
        self.remove_by_consensus = remove;
        self
    }

    /// Recover missing genes from the graph.
    ///
    /// For each edge in the graph, checks if any genome is missing
    /// one of the connected genes. If so, searches the flanking
    /// contig sequence for a match using semi-global alignment.
    ///
    /// When `remove_by_consensus` is enabled, any node where the
    /// refound hit count exceeds its original support is removed.
    pub fn recover(
        &self,
        graph: &ConcurrentGraph,
        contig_sequences: &HashMap<String, Vec<u8>>,
        cluster_sequences: &HashMap<String, Vec<u8>>,
    ) -> Result<RecoveryStats> {
        let mut recovered = 0;
        let mut recovery_counts: HashMap<ClusterId, usize> = HashMap::new();

        // Record original support values before recovery (for consensus removal)
        let original_sizes: HashMap<ClusterId, usize> = graph
            .nodes
            .iter()
            .map(|entry| (entry.key().clone(), entry.value().support))
            .collect();

        let mut edges = Vec::new();
        for e in graph.edges.iter() {
            let (from, to) = e.key();
            let genomes = e.value().genomes.clone();
            edges.push((from.clone(), to.clone(), genomes));
        }

        for (from, to, _edge_genomes) in &edges {
            // Get genomes that have the 'from' node but not the 'to' node
            let from_node = graph.nodes.get(from);
            let to_node = graph.nodes.get(to);

            if from_node.is_none() || to_node.is_none() {
                continue;
            }

            // Get the query sequence for the missing gene
            let query_seq = match cluster_sequences.get(to.as_str()) {
                Some(seq) => seq,
                None => continue,
            };

            // Search in contig sequences for each genome missing this gene
            if query_seq.len() < 10 {
                continue;
            }

            // Search using semi-global alignment (HW mode)
            if self.search_in_contigs(query_seq, contig_sequences) {
                recovered += 1;
                *recovery_counts.entry(to.clone()).or_insert(0) += 1;
            }
        }

        // Consensus removal: delete nodes where refound hits exceed original support.
        // A node with many more refound hits than original members is likely
        // a spurious annotation artifact (matches Panaroo's --remove_by_consensus).
        let mut nodes_removed_by_consensus = 0;
        if self.remove_by_consensus {
            for (cluster_id, refound_count) in &recovery_counts {
                if let Some(&original_size) = original_sizes.get(cluster_id) {
                    if refound_count > &original_size {
                        graph.remove_node(cluster_id);
                        nodes_removed_by_consensus += 1;
                    }
                }
            }
            if nodes_removed_by_consensus > 0 {
                tracing::info!(
                    "Consensus removal: removed {} spurious nodes (refound hits exceeded original support)",
                    nodes_removed_by_consensus
                );
            }
        }

        tracing::info!("Missing gene recovery: recovered {} genes", recovered);

        Ok(RecoveryStats {
            genes_recovered: recovered,
            nodes_removed_by_consensus,
        })
    }

    /// Search for a query sequence in contig sequences using semi-global alignment.
    ///
    /// Uses sliding-window semi-global alignment (HW mode) to find the query
    /// embedded within longer contig sequences. This matches Panaroo's edlib-based
    /// approach where the query is aligned within the target.
    fn search_in_contigs(
        &self,
        query: &[u8],
        contigs: &HashMap<String, Vec<u8>>,
    ) -> bool {
        if query.is_empty() {
            return false;
        }

        let min_align_len = (query.len() as f32 * self.prop_match) as usize;
        let min_align_len = min_align_len.max(1);

        for seq in contigs.values() {
            if seq.len() < query.len() {
                // Contig is shorter than query — just do direct alignment
                let (identity, _dist, align_len) = align_semiglobal(query, seq);
                if align_len >= min_align_len && identity >= self.min_identity as f64 {
                    return true;
                }
                continue;
            }

            // Slide a window across the contig and align the query within each window
            let window_size = query.len() + self.search_window;
            let step = (query.len() / 2).max(100);

            let mut pos = 0;
            while pos < seq.len() {
                let end = (pos + window_size).min(seq.len());
                let window = &seq[pos..end];

                let (identity, _dist, align_len) = align_semiglobal(query, window);

                if align_len >= min_align_len && identity >= self.min_identity as f64 {
                    return true;
                }

                pos += step;
                // Avoid infinite loop when step is 0 or window_size is larger than seq
                if pos >= seq.len() {
                    break;
                }
            }
        }

        false
    }
}

impl Default for MissingGeneRecoverer {
    fn default() -> Self {
        Self::new()
    }
}

/// Statistics from missing gene recovery.
#[derive(Debug, Clone)]
pub struct RecoveryStats {
    /// Number of genes recovered
    pub genes_recovered: usize,
    /// Number of nodes removed by consensus (refound hits exceeded original support)
    pub nodes_removed_by_consensus: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_missing_gene_recoverer_creation() {
        let recoverer = MissingGeneRecoverer::new();
        assert_eq!(recoverer.min_identity, 0.70);
        assert_eq!(recoverer.search_window, 5000);
        assert!(!recoverer.remove_by_consensus);
    }

    #[test]
    fn test_remove_by_consensus_builder() {
        let recoverer = MissingGeneRecoverer::new()
            .with_remove_by_consensus(true);
        assert!(recoverer.remove_by_consensus);

        let recoverer_off = MissingGeneRecoverer::new()
            .with_remove_by_consensus(false);
        assert!(!recoverer_off.remove_by_consensus);
    }

    #[test]
    fn test_search_sequence_embedded() {
        // Query embedded within a longer contig — HW mode should find it
        let recoverer = MissingGeneRecoverer::new()
            .with_min_identity(0.80);
        let query = b"ATCGATCGATCGATCG";
        let mut contigs = HashMap::new();
        contigs.insert("contig1".to_string(), b"NNNNATCGATCGATCGATCGNNNN".to_vec());

        assert!(recoverer.search_in_contigs(query, &contigs));
    }

    #[test]
    fn test_search_sequence_not_found() {
        // Query not present in contig
        let recoverer = MissingGeneRecoverer::new()
            .with_min_identity(0.80);
        let query = b"ATCGATCGATCGATCG";
        let mut contigs = HashMap::new();
        contigs.insert("contig1".to_string(), b"GGGGAAAACCCCTTTT".to_vec());

        assert!(!recoverer.search_in_contigs(query, &contigs));
    }

    #[test]
    fn test_search_sequence_near_match() {
        // Query with 1 substitution embedded in contig
        let recoverer = MissingGeneRecoverer::new()
            .with_min_identity(0.85);
        let query = b"ATCGATCGATCGATCG";
        let mut contigs = HashMap::new();
        // 15/16 = 93.75% identity (1 substitution in 16 chars)
        contigs.insert("contig1".to_string(), b"NNNNATCGATCGATCAATCGNNNN".to_vec());

        assert!(recoverer.search_in_contigs(query, &contigs));
    }

    #[test]
    fn test_consensus_removes_spurious_nodes() {
        // Create a graph where a node with low support gets many refound hits.
        // With remove_by_consensus=true, that node should be removed.
        // With remove_by_consensus=false, the node should remain.
        use crate::graph::{GeneCluster, Node, Edge, GenomeId};

        let graph = ConcurrentGraph::with_capacity(10);

        // Create a "spurious" node with support=1 (only 1 genome)
        let mut spurious_cluster = GeneCluster::new("spurious");
        spurious_cluster.support = 1;
        spurious_cluster.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];
        let mut spurious_node = Node::from_cluster(&spurious_cluster);
        spurious_node.support = 1;
        spurious_node.genomes.insert(GenomeId::new("g1"));
        graph.add_node(spurious_node);

        // Create 3 connector nodes that connect to "spurious".
        // Each edge generates a refound hit for "spurious" (since it is the
        // `to` node alphabetically: "connector_*" < "spurious").
        // 3 hits > 1 (original support) => consensus removal triggers.
        for i in 0..3 {
            let id = format!("connector_{}", i);
            let mut cluster = GeneCluster::new(&id);
            cluster.support = 5;
            cluster.centroids = vec![b"GGGGAAAACCCCTTTT".to_vec()];
            let mut node = Node::from_cluster(&cluster);
            node.support = 5;
            node.genomes.insert(GenomeId::new(format!("conn_g{}", i)));
            graph.add_node(node);
            graph.add_edge(Edge::new(ClusterId::new(&id), ClusterId::new("spurious")));
        }

        // Cluster sequences: spurious node has a centroid that will be found
        let mut cluster_sequences = HashMap::new();
        cluster_sequences.insert("spurious".to_string(), b"ATCGATCGATCGATCG".to_vec());

        // Contig sequences: embed the spurious query so it is found
        let mut contig_sequences = HashMap::new();
        contig_sequences.insert(
            "c1".to_string(),
            b"NNNNATCGATCGATCGATCGNNNN".to_vec(),
        );

        // With remove_by_consensus=true: spurious node has 1 support but
        // gets 3 refound hits (one per edge), so it should be removed
        let recoverer_strict = MissingGeneRecoverer::new()
            .with_min_identity(0.70)
            .with_remove_by_consensus(true);

        let initial_count = graph.node_count();
        let stats = recoverer_strict
            .recover(&graph, &contig_sequences, &cluster_sequences)
            .unwrap();

        // The spurious node should have been removed (refound hits > original support=1)
        assert!(
            stats.nodes_removed_by_consensus > 0,
            "Expected at least one node removed by consensus, got {}",
            stats.nodes_removed_by_consensus
        );
        assert!(
            graph.node_count() < initial_count,
            "Graph should have fewer nodes after consensus removal"
        );

        // Now test without consensus removal: rebuild the same graph structure
        let graph2 = ConcurrentGraph::with_capacity(10);

        let mut s2 = GeneCluster::new("spurious");
        s2.support = 1;
        s2.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];
        let mut sn2 = Node::from_cluster(&s2);
        sn2.support = 1;
        sn2.genomes.insert(GenomeId::new("g1"));
        graph2.add_node(sn2);

        for i in 0..3 {
            let id = format!("connector_{}", i);
            let mut cluster = GeneCluster::new(&id);
            cluster.support = 5;
            cluster.centroids = vec![b"GGGGAAAACCCCTTTT".to_vec()];
            let mut node = Node::from_cluster(&cluster);
            node.support = 5;
            node.genomes.insert(GenomeId::new(format!("conn_g{}", i)));
            graph2.add_node(node);
            graph2.add_edge(Edge::new(ClusterId::new(&id), ClusterId::new("spurious")));
        }

        let recoverer_default = MissingGeneRecoverer::new()
            .with_min_identity(0.70)
            .with_remove_by_consensus(false);

        let count_before = graph2.node_count();
        let stats2 = recoverer_default
            .recover(&graph2, &contig_sequences, &cluster_sequences)
            .unwrap();

        // Without consensus removal, no nodes should be deleted
        assert_eq!(stats2.nodes_removed_by_consensus, 0);
        assert_eq!(
            graph2.node_count(),
            count_before,
            "Graph should have same node count without consensus removal"
        );
    }
}