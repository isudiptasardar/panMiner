//! Misassembly edge cleaning.
//!
//! Removes edges likely caused by misassemblies — low-support edges
//! near contig ends or disproportionately weak edges.

use crate::error::Result;
use crate::graph::ConcurrentGraph;
use crate::graph::ClusterId;

/// Cleans misassembly edges from the pangenome graph.
pub struct MisassemblyEdgeCleaner {
    /// Minimum edge support threshold
    edge_support_threshold: usize,
    /// Proportion threshold for relative weakness (hardcoded at 0.05 like Panaroo)
    proportion_threshold: f64,
}

/// Statistics from edge cleaning.
#[derive(Debug, Clone)]
pub struct CleaningStats {
    /// Number of edges removed
    pub edges_removed: usize,
}

impl MisassemblyEdgeCleaner {
    /// Create a new cleaner with the given threshold.
    pub fn new(edge_support_threshold: usize) -> Self {
        Self {
            edge_support_threshold,
            proportion_threshold: 0.05,
        }
    }

    /// Create from correction mode and genome count.
    ///
    /// Thresholds follow Panaroo conventions:
    /// - Strict/Default: max(2, ceil(1% × n))
    /// - Sensitive: 0 (disabled)
    pub fn from_mode(mode: &crate::config::CorrectionMode, num_genomes: usize) -> Self {
        let threshold = match mode {
            crate::config::CorrectionMode::Strict => {
                std::cmp::max(2, (num_genomes as f64 * 0.01).ceil() as usize)
            }
            crate::config::CorrectionMode::Default => {
                std::cmp::max(2, (num_genomes as f64 * 0.01).ceil() as usize)
            }
            crate::config::CorrectionMode::Sensitive => 0, // disabled
        };
        Self::new(threshold)
    }

    /// Clean misassembly edges from the graph.
    ///
    /// Two removal criteria (matching Panaroo's clean_misassembly_edges):
    /// 1. Edges connected to contig-end nodes with support < threshold
    /// 2. Edges with support < 5% of the smaller node's support AND < threshold
    pub fn clean(&self, graph: &ConcurrentGraph) -> Result<CleaningStats> {
        if self.edge_support_threshold == 0 {
            return Ok(CleaningStats { edges_removed: 0 });
        }

        let mut bad_edges: Vec<(ClusterId, ClusterId)> = Vec::new();

        // Criterion 1: Edges near contig-end nodes with low support
        for entry in graph.nodes.iter() {
            let node = entry.value();
            if !node.contig_end_genomes.is_empty() {
                let cluster_id = entry.key();
                // Use adjacency index for O(degree) lookup instead of O(E) scan
                for neighbor in graph.neighbors(cluster_id) {
                    let key = if cluster_id < &neighbor {
                        (cluster_id.clone(), neighbor.clone())
                    } else {
                        (neighbor.clone(), cluster_id.clone())
                    };
                    if let Some(edge) = graph.edges.get(&key) {
                        if edge.support < self.edge_support_threshold {
                            bad_edges.push(key);
                        }
                    }
                }
            }
        }

        // Criterion 2: Disproportionately weak edges
        for edge_entry in graph.edges.iter() {
            let (from, to) = edge_entry.key();
            let edge = edge_entry.value();
            let support_a = graph.nodes.get(from).map(|n| n.support).unwrap_or(0);
            let support_b = graph.nodes.get(to).map(|n| n.support).unwrap_or(0);
            let min_node_support = support_a.min(support_b);

            if min_node_support > 0 {
                let proportion = edge.support as f64 / min_node_support as f64;
                if proportion < self.proportion_threshold
                    && edge.support < self.edge_support_threshold
                {
                    bad_edges.push((from.clone(), to.clone()));
                }
            }
        }

        // Deduplicate
        bad_edges.sort();
        bad_edges.dedup();

        // Remove bad edges (using remove_edge to keep adjacency index consistent)
        let edges_removed = bad_edges.len();
        for (from, to) in &bad_edges {
            graph.remove_edge(from, to);
        }

        tracing::info!(
            "Misassembly edge cleaning: removed {} edges",
            edges_removed
        );

        Ok(CleaningStats { edges_removed })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Node, GeneCluster, GenomeId};

    #[test]
    fn test_misassembly_edge_cleaning_contig_end() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Add contig-end node with low support
        let mut end_node = Node::from_cluster(&{
            let mut c = GeneCluster::new("end");
            c.support = 1;
            c
        });
        end_node.contig_end_genomes.insert(GenomeId::new("g1"));
        graph.add_node(end_node);

        // Add high-support core node
        let core_node = Node::from_cluster(&{
            let mut c = GeneCluster::new("core");
            c.support = 100;
            c
        });
        graph.add_node(core_node);

        // Add edge with low support
        graph.add_edge_genome(
            ClusterId::new("end"),
            ClusterId::new("core"),
            GenomeId::new("g1"),
        );

        let cleaner = MisassemblyEdgeCleaner::new(2);
        let stats = cleaner.clean(&graph).unwrap();

        assert_eq!(stats.edges_removed, 1);
        assert_eq!(graph.edge_count(), 0);
    }

    #[test]
    fn test_misassembly_edge_cleaning_disproportionate() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Two high-support nodes
        let node_a = Node::from_cluster(&{
            let mut c = GeneCluster::new("a");
            c.support = 100;
            c
        });
        let node_b = Node::from_cluster(&{
            let mut c = GeneCluster::new("b");
            c.support = 80;
            c
        });
        graph.add_node(node_a);
        graph.add_node(node_b);

        // Edge with support 2 (2.5% of min node support 80)
        graph.add_edge_genome(ClusterId::new("a"), ClusterId::new("b"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("a"), ClusterId::new("b"), GenomeId::new("g2"));

        let cleaner = MisassemblyEdgeCleaner::new(10);
        let stats = cleaner.clean(&graph).unwrap();

        // Edge support (2) < 5% of 80 AND < 10 → should be removed
        assert_eq!(stats.edges_removed, 1);
    }

    #[test]
    fn test_misassembly_sensitive_disabled() {
        let graph = ConcurrentGraph::with_capacity(10);
        let cleaner = MisassemblyEdgeCleaner::from_mode(
            &crate::config::CorrectionMode::Sensitive, 100
        );
        let stats = cleaner.clean(&graph).unwrap();
        assert_eq!(stats.edges_removed, 0);
    }

    #[test]
    fn test_misassembly_preserves_strong_edges() {
        let graph = ConcurrentGraph::with_capacity(10);

        let node_a = Node::from_cluster(&{
            let mut c = GeneCluster::new("a");
            c.support = 50;
            c
        });
        let node_b = Node::from_cluster(&{
            let mut c = GeneCluster::new("b");
            c.support = 50;
            c
        });
        graph.add_node(node_a);
        graph.add_node(node_b);

        // Strong edge: 30/50 = 60% — should NOT be removed
        for i in 0..30 {
            graph.add_edge_genome(
                ClusterId::new("a"),
                ClusterId::new("b"),
                GenomeId::new(format!("g{}", i)),
            );
        }

        let cleaner = MisassemblyEdgeCleaner::new(2);
        let stats = cleaner.clean(&graph).unwrap();

        assert_eq!(stats.edges_removed, 0);
        assert_eq!(graph.edge_count(), 1);
    }
}