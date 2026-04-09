//! Contamination removal from the pangenome graph.
//!
//! Recursively removes low-support nodes with degree <= 1,
//! which typically represent contamination or spurious annotations.

use crate::config::CorrectionMode;
use crate::error::Result;
use crate::graph::ConcurrentGraph;

/// Removes contamination from the pangenome graph.
///
/// Low-support nodes with degree <= 1 are recursively removed.
/// The threshold depends on the correction mode.
pub struct ContaminationRemover {
    /// Support threshold below which nodes are removed
    threshold: usize,
    /// Maximum iterations for recursive removal
    max_iterations: usize,
}

impl ContaminationRemover {
    /// Create a new contamination remover.
    pub fn new(threshold: usize) -> Self {
        Self {
            threshold,
            max_iterations: 100,
        }
    }

    /// Create from correction mode.
    pub fn from_mode(mode: &CorrectionMode, num_genomes: usize) -> Self {
        let threshold = match mode {
            CorrectionMode::Strict => (num_genomes as f64 * 0.05).ceil() as usize,
            CorrectionMode::Default => 2,
            CorrectionMode::Sensitive => 1,
        };

        Self::new(threshold)
    }

    /// Remove contamination from the graph.
    ///
    /// Recursively removes low-support degree-1 nodes until
    /// no more can be removed.
    pub fn remove(&self, graph: &ConcurrentGraph) -> Result<RemovalStats> {
        let mut total_removed = 0;
        let mut iteration = 0;

        loop {
            if iteration >= self.max_iterations {
                break;
            }

            let to_remove = graph.find_low_support_nodes(self.threshold);

            if to_remove.is_empty() {
                break;
            }

            let removed_count = to_remove.len();
            graph.remove_nodes_parallel(&to_remove);
            total_removed += removed_count;
            iteration += 1;

            tracing::debug!(
                "Contamination removal iteration {}: removed {} nodes",
                iteration,
                removed_count
            );
        }

        tracing::info!(
            "Contamination removal: removed {} nodes in {} iterations",
            total_removed,
            iteration
        );

        Ok(RemovalStats {
            nodes_removed: total_removed,
            iterations: iteration,
        })
    }
}

/// Statistics from contamination removal.
#[derive(Debug, Clone)]
pub struct RemovalStats {
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
    fn test_contamination_removal() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Add a low-support degree-1 node
        let low = Node::from_cluster(&{
            let mut c = GeneCluster::new("low");
            c.support = 1;
            c
        });
        graph.add_node(low);

        // Add a high-support node
        let high = Node::from_cluster(&{
            let mut c = GeneCluster::new("high");
            c.support = 100;
            c
        });
        graph.add_node(high);

        // Connect them
        graph.add_edge_genome(
            ClusterId::new("low"),
            ClusterId::new("high"),
            GenomeId::new("g1"),
        );

        let remover = ContaminationRemover::new(2);
        let stats = remover.remove(&graph).unwrap();

        assert_eq!(stats.nodes_removed, 1);
        assert_eq!(graph.node_count(), 1);
    }
}
