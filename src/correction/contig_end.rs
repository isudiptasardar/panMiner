//! Contig-end pruning module for PanMiner.
//!
//! Prunes nodes that represent contig ends (terminal nodes with degree 1)
//! that have support below a minimum threshold.

use crate::config::CorrectionMode;
use crate::error::Result;
use crate::graph::ConcurrentGraph;

/// Prunes contig-end nodes from the pangenome graph.
///
/// Contig-end nodes are terminal nodes (degree 1) that likely represent
/// incomplete genes at contig boundaries. Low-support contig-end nodes
/// are removed to improve graph quality.
pub struct ContigEndPruner {
    /// Minimum support threshold for keeping a contig-end node
    min_support: usize,
    /// Maximum iterations for recursive pruning
    max_iterations: usize,
}

impl ContigEndPruner {
    /// Create a new contig-end pruner with default settings.
    pub fn new() -> Self {
        Self {
            min_support: 1,
            max_iterations: 100,
        }
    }

    /// Create from correction mode and genome count.
    ///
    /// Thresholds follow Panaroo conventions:
    /// - Strict: max(2, ceil(0.05 × n))
    /// - Default: max(2, ceil(0.01 × n))
    /// - Sensitive: 2 (minimal removal)
    pub fn from_mode(mode: &CorrectionMode, num_genomes: usize) -> Self {
        let min_support = match mode {
            CorrectionMode::Strict => std::cmp::max(2, (num_genomes as f64 * 0.05).ceil() as usize),
            CorrectionMode::Default => std::cmp::max(2, (num_genomes as f64 * 0.01).ceil() as usize),
            CorrectionMode::Sensitive => 2,
        };
        Self {
            min_support,
            max_iterations: 100,
        }
    }

    /// Set the minimum support threshold.
    pub fn with_min_support(mut self, min_support: usize) -> Self {
        self.min_support = min_support;
        self
    }

    /// Prune contig-end nodes from the graph.
    ///
    /// Removes nodes that:
    /// - Have degree 1 (only one connected edge)
    /// - Have support below the minimum threshold
    pub fn prune(&self, graph: &ConcurrentGraph) -> Result<PruningStats> {
        let mut nodes_removed = 0;
        let mut iteration = 0;

        loop {
            if iteration >= self.max_iterations {
                break;
            }

            let to_remove = graph
                .nodes
                .iter()
                .filter(|entry| {
                    let node = entry.value();
                    // Check if it's a contig-end node (degree 1) with low support
                    !node.contig_end_genomes.is_empty() && graph.is_degree_one(entry.key()) && node.support < self.min_support
                })
                .map(|entry| entry.key().clone())
                .collect::<Vec<_>>();

            if to_remove.is_empty() {
                break;
            }

            let removed_count = to_remove.len();
            graph.remove_nodes_parallel(&to_remove);
            nodes_removed += removed_count;
            iteration += 1;

            tracing::debug!(
                "Contig-end pruning iteration {}: removed {} nodes",
                iteration,
                removed_count
            );
        }

        tracing::info!(
            "Contig-end pruning: removed {} nodes in {} iterations",
            nodes_removed,
            iteration
        );

        Ok(PruningStats {
            nodes_removed,
            iterations: iteration,
        })
    }
}

impl Default for ContigEndPruner {
    fn default() -> Self {
        Self::new()
    }
}

/// Statistics from contig-end pruning.
#[derive(Debug, Clone)]
pub struct PruningStats {
    /// Total nodes removed
    pub nodes_removed: usize,
    /// Number of iterations
    pub iterations: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Node, GeneCluster, GenomeId, ClusterId};

    #[test]
    fn test_contig_end_pruning() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Add a low-support contig-end node
        let end_node = Node::from_cluster(&{
            let mut c = GeneCluster::new("end");
            c.support = 1;
            c
        });
        let mut end_node = end_node;
        end_node.contig_end_genomes.insert(GenomeId::new("g1"));
        graph.add_node(end_node);

        // Add a high-support node
        let core_node = Node::from_cluster(&{
            let mut c = GeneCluster::new("core");
            c.support = 100;
            c
        });
        graph.add_node(core_node);

        // Connect them (makes end node degree-1)
        graph.add_edge_genome(
            ClusterId::new("end"),
            ClusterId::new("core"),
            GenomeId::new("g1"),
        );

        let pruner = ContigEndPruner::new().with_min_support(2);
        let stats = pruner.prune(&graph).unwrap();

        assert_eq!(stats.nodes_removed, 1);
        assert_eq!(graph.node_count(), 1);
    }

    #[test]
    fn test_contig_end_not_pruned_if_high_support() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Add a high-support contig-end node
        let end_node = Node::from_cluster(&{
            let mut c = GeneCluster::new("end");
            c.support = 5;
            c
        });
        let mut end_node = end_node;
        end_node.contig_end_genomes.insert(GenomeId::new("g1"));
        graph.add_node(end_node);

        // Add another node
        let other_node = Node::from_cluster(&{
            let mut c = GeneCluster::new("other");
            c.support = 100;
            c
        });
        graph.add_node(other_node);

        // Connect them
        graph.add_edge_genome(
            ClusterId::new("end"),
            ClusterId::new("other"),
            GenomeId::new("g1"),
        );

        let pruner = ContigEndPruner::new().with_min_support(2);
        let stats = pruner.prune(&graph).unwrap();

        // Should not remove because support >= threshold
        assert_eq!(stats.nodes_removed, 0);
        assert_eq!(graph.node_count(), 2);
    }

    #[test]
    fn test_contig_end_not_pruned_if_not_degree_one() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Add a contig-end node connected to multiple nodes
        let end_node = Node::from_cluster(&{
            let mut c = GeneCluster::new("end");
            c.support = 1;
            c
        });
        let mut end_node = end_node;
        end_node.contig_end_genomes.insert(GenomeId::new("g1"));
        graph.add_node(end_node);

        // Add two other nodes
        let node1 = Node::from_cluster(&{
            let mut c = GeneCluster::new("node1");
            c.support = 100;
            c
        });
        graph.add_node(node1);

        let node2 = Node::from_cluster(&{
            let mut c = GeneCluster::new("node2");
            c.support = 100;
            c
        });
        graph.add_node(node2);

        // Connect end node to both (degree 2)
        graph.add_edge_genome(
            ClusterId::new("end"),
            ClusterId::new("node1"),
            GenomeId::new("g1"),
        );
        graph.add_edge_genome(
            ClusterId::new("end"),
            ClusterId::new("node2"),
            GenomeId::new("g1"),
        );

        let pruner = ContigEndPruner::new().with_min_support(2);
        let stats = pruner.prune(&graph).unwrap();

        // Should not remove because degree > 1
        assert_eq!(stats.nodes_removed, 0);
        assert_eq!(graph.node_count(), 3);
    }

    #[test]
    fn test_from_mode_thresholds() {
        // Strict: max(2, ceil(0.05 * n))
        let strict = ContigEndPruner::from_mode(&CorrectionMode::Strict, 100);
        assert_eq!(strict.min_support, 5); // ceil(0.05 * 100) = 5

        let strict_small = ContigEndPruner::from_mode(&CorrectionMode::Strict, 10);
        assert_eq!(strict_small.min_support, 2); // max(2, ceil(0.5)) = 2

        // Default: max(2, ceil(0.01 * n))
        let default_mode = ContigEndPruner::from_mode(&CorrectionMode::Default, 100);
        assert_eq!(default_mode.min_support, 2); // max(2, ceil(1)) = 2

        let default_mode_large = ContigEndPruner::from_mode(&CorrectionMode::Default, 500);
        assert_eq!(default_mode_large.min_support, 5); // max(2, ceil(5)) = 5

        // Sensitive: always 2
        let sensitive = ContigEndPruner::from_mode(&CorrectionMode::Sensitive, 100);
        assert_eq!(sensitive.min_support, 2);

        let sensitive_large = ContigEndPruner::from_mode(&CorrectionMode::Sensitive, 1000);
        assert_eq!(sensitive_large.min_support, 2);
    }
}
