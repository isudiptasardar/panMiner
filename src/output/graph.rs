//! GML graph output for Cytoscape visualization.

use std::io::Write;
use std::path::Path;

use crate::error::Result;
use crate::graph::PangenomeGraph;

/// Writer for GML (Graph Modelling Language) format.
///
/// GML files can be visualized with Cytoscape and other graph tools.
pub struct GmlWriter;

impl GmlWriter {
    /// Write the pangenome graph to GML format.
    pub fn write(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        writeln!(file, "graph [")?;
        writeln!(file, "  directed 0")?;

        // Write nodes
        for (cluster_id, node) in &graph.nodes {
            writeln!(file, "  node [")?;
            writeln!(file, "    id \"{}\"", cluster_id)?;
            writeln!(file, "    label \"{}\"", cluster_id)?;
            writeln!(file, "    support {}", node.support)?;
            writeln!(file, "    is_paralog {}", if node.is_paralog { 1 } else { 0 })?;

            if let Some(annotation) = node.annotations.iter().next() {
                // Escape quotes in annotation
                let escaped = annotation.replace('"', "\\\"");
                writeln!(file, "    annotation \"{}\"", escaped)?;
            }

            writeln!(file, "  ]")?;
        }

        // Write edges
        for ((from, to), edge) in &graph.edges {
            writeln!(file, "  edge [")?;
            writeln!(file, "    source \"{}\"", from)?;
            writeln!(file, "    target \"{}\"", to)?;
            writeln!(file, "    support {}", edge.support)?;
            writeln!(file, "    genomes {}", edge.genomes.len())?;
            writeln!(file, "  ]")?;
        }

        writeln!(file, "]")?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Node, Edge, GeneCluster, ClusterId, GenomeId};
    use tempfile::NamedTempFile;

    #[test]
    fn test_gml_output() {
        let mut graph = PangenomeGraph::new();

        let node = Node::from_cluster(&{
            let mut c = GeneCluster::new("c1");
            c.support = 5;
            c
        });
        graph.add_node(node);

        let node2 = Node::from_cluster(&{
            let mut c = GeneCluster::new("c2");
            c.support = 3;
            c
        });
        graph.add_node(node2);

        let mut edge = Edge::new(ClusterId::new("c1"), ClusterId::new("c2"));
        edge.add_genome(GenomeId::new("g1"));
        graph.add_edge(edge);

        let temp = NamedTempFile::new().unwrap();
        GmlWriter::write(&graph, temp.path()).unwrap();

        let content = std::fs::read_to_string(temp.path()).unwrap();
        assert!(content.contains("graph ["));
        assert!(content.contains("node ["));
        assert!(content.contains("edge ["));
    }
}
