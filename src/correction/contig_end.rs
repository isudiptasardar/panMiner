//! Contig-end pruning module for PanMiner.
//!
//! Prunes nodes that represent contig ends (terminal nodes with degree 1)
//! that have support below a minimum threshold.

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
}

impl ContigEndPruner {
    /// Create a new contig-end pruner with default settings.
    pub fn new() -> Self {
        Self {
            min_support: 1,
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
}
