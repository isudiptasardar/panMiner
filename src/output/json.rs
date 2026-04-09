//! JSON/JSONL output for programmatic access.

use std::io::Write;
use std::path::Path;
use serde::Serialize;

use crate::error::Result;
use crate::graph::{PangenomeGraph, BitPackedMatrix};

/// Writer for JSON output.
pub struct JsonWriter;

/// JSON representation of the pangenome summary.
#[derive(Serialize)]
struct PangenomeSummary {
    version: String,
    num_genomes: usize,
    num_clusters: usize,
    num_core: usize,
    num_accessory: usize,
    num_edges: usize,
    clusters: Vec<ClusterSummary>,
}

/// JSON representation of a cluster.
#[derive(Serialize)]
struct ClusterSummary {
    id: String,
    support: usize,
    is_paralog: bool,
    annotations: Vec<String>,
    neighbors: Vec<String>,
}

impl JsonWriter {
    /// Write the pangenome graph and matrix to JSON.
    pub fn write(
        graph: &PangenomeGraph,
        _matrix: &BitPackedMatrix,
        path: &Path,
    ) -> Result<()> {
        let total_genomes = graph.genomes.len().max(1);
        let core_threshold = (total_genomes as f32 * 0.99).ceil() as usize;

        let mut clusters: Vec<ClusterSummary> = Vec::new();
        let mut num_core = 0;

        for (cluster_id, node) in &graph.nodes {
            if node.support >= core_threshold {
                num_core += 1;
            }

            let neighbors: Vec<String> = graph.neighbors(cluster_id)
                .iter()
                .map(|id| id.to_string())
                .collect();

            clusters.push(ClusterSummary {
                id: cluster_id.to_string(),
                support: node.support,
                is_paralog: node.is_paralog,
                annotations: node.annotations.iter().cloned().collect(),
                neighbors,
            });
        }

        // Sort clusters by support (descending)
        clusters.sort_by(|a, b| b.support.cmp(&a.support));

        let summary = PangenomeSummary {
            version: crate::VERSION.to_string(),
            num_genomes: total_genomes,
            num_clusters: graph.node_count(),
            num_core,
            num_accessory: graph.node_count() - num_core,
            num_edges: graph.edge_count(),
            clusters,
        };

        let json = serde_json::to_string_pretty(&summary)?;
        std::fs::write(path, json)?;

        Ok(())
    }

    /// Write the pangenome to JSONL (one JSON object per line).
    ///
    /// JSONL is more suitable for streaming and large datasets.
    pub fn write_jsonl(
        graph: &PangenomeGraph,
        path: &Path,
    ) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        for (cluster_id, node) in &graph.nodes {
            let neighbors: Vec<String> = graph.neighbors(cluster_id)
                .iter()
                .map(|id| id.to_string())
                .collect();

            let entry = ClusterSummary {
                id: cluster_id.to_string(),
                support: node.support,
                is_paralog: node.is_paralog,
                annotations: node.annotations.iter().cloned().collect(),
                neighbors,
            };

            let line = serde_json::to_string(&entry)?;
            writeln!(file, "{}", line)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Node, GeneCluster};
    use tempfile::NamedTempFile;

    #[test]
    fn test_json_output() {
        let mut graph = PangenomeGraph::new();

        let node = Node::from_cluster(&{
            let mut c = GeneCluster::new("c1");
            c.support = 5;
            c
        });
        graph.add_node(node);

        let matrix = BitPackedMatrix::new(5, 1);

        let temp = NamedTempFile::new().unwrap();
        JsonWriter::write(&graph, &matrix, temp.path()).unwrap();

        let content = std::fs::read_to_string(temp.path()).unwrap();
        assert!(content.contains("num_clusters"));
        assert!(content.contains("c1"));
    }
}
