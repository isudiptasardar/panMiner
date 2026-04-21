//! GML graph output for Cytoscape visualization.

use std::io::Write;
use std::path::Path;

use crate::error::Result;
use crate::graph::PangenomeGraph;

/// Escape a string for GML format (handle quotes and backslashes).
fn escape_gml_string(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

/// Writer for GML (Graph Modelling Language) format.
///
/// GML files can be visualized with Cytoscape and other graph tools.
pub struct GmlWriter;

impl GmlWriter {
    /// Write the pangenome graph to GML format.
    pub fn write(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;
        let mut writer = std::io::BufWriter::new(&mut file);

        writeln!(writer, "graph [")?;
        writeln!(writer, "  directed 0")?;

        // Nodes with full attributes (Panaroo-compatible)
        for (_id, node) in &graph.nodes {
            writeln!(writer, "  node [")?;
            writeln!(writer, "    id \"{}\"", node.cluster_id)?;
            writeln!(writer, "    label \"{}\"", node.cluster_id)?;
            writeln!(writer, "    support {}", node.support)?;
            writeln!(writer, "    is_paralog {}", if node.is_paralog { 1 } else { 0 })?;
            writeln!(writer, "    is_highly_variable {}", if node.is_highly_variable { 1 } else { 0 })?;

            // Length of centroid sequence
            let length = node.centroid_sequences.first().map(|s| s.len()).unwrap_or(0);
            writeln!(writer, "    length {}", length)?;

            // Centroid DNA sequences (comma-separated for multi-centroid nodes)
            if !node.centroid_sequences.is_empty() {
                let centroid_seqs: Vec<String> = node.centroid_sequences.iter()
                    .map(|s| String::from_utf8_lossy(s).to_string())
                    .collect();
                writeln!(writer, "    seq \"{}\"", escape_gml_string(&centroid_seqs.join(",")))?;

                // Protein sequence from first centroid
                if let Some(seq) = node.centroid_sequences.first() {
                    let protein = crate::io::translate(seq);
                    let protein_str = String::from_utf8_lossy(&protein);
                    if !protein_str.is_empty() {
                        writeln!(writer, "    protein \"{}\"", escape_gml_string(&protein_str))?;
                    }
                }
            }

            // Genome IDs (bracket list format for GML compatibility)
            let genome_ids: Vec<String> = node.genomes.iter()
                .map(|g| g.as_str().to_string())
                .collect();
            if !genome_ids.is_empty() {
                writeln!(writer, "    genomes [")?;
                for gid in &genome_ids {
                    writeln!(writer, "      \"{}\"", gid)?;
                }
                writeln!(writer, "    ]")?;
            }

            // Gene members (semicolon-separated)
            let all_members: Vec<String> = node.gene_members.values()
                .flatten()
                .cloned()
                .collect();
            if !all_members.is_empty() {
                writeln!(writer, "    member \"{}\"", all_members.join(";"))?;
            }

            // Contig-end genomes (comma-separated)
            let contig_end_ids: Vec<String> = node.contig_end_genomes.iter()
                .map(|g| g.as_str().to_string())
                .collect();
            if !contig_end_ids.is_empty() {
                writeln!(writer, "    contig_end_genomes \"{}\"", contig_end_ids.join(","))?;
            }

            // Annotation
            if let Some(ann) = node.annotations.iter().next() {
                writeln!(writer, "    annotation \"{}\"", escape_gml_string(ann))?;
            }

            writeln!(writer, "  ]")?;
        }

        // Edges with genome IDs
        for (_key, edge) in &graph.edges {
            writeln!(writer, "  edge [")?;
            writeln!(writer, "    source \"{}\"", edge.from)?;
            writeln!(writer, "    target \"{}\"", edge.to)?;
            writeln!(writer, "    support {}", edge.support)?;

            // Genome IDs on edges (bracket list format for GML compatibility)
            let edge_genome_ids: Vec<String> = edge.genomes.iter()
                .map(|g| g.as_str().to_string())
                .collect();
            if !edge_genome_ids.is_empty() {
                writeln!(writer, "    genomes [")?;
                for gid in &edge_genome_ids {
                    writeln!(writer, "      \"{}\"", gid)?;
                }
                writeln!(writer, "    ]")?;
            }

            writeln!(writer, "  ]")?;
        }

        writeln!(writer, "]")?;
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

    #[test]
    fn test_gml_output_with_full_attributes() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.gml");

        let mut graph = PangenomeGraph::new();
        let mut node = Node::from_cluster(&GeneCluster::new("c1"));
        node.centroid_sequences = vec![b"ATGCGT".to_vec()];
        node.support = 3;
        node.genomes.insert(GenomeId::new("genome1"));
        node.genomes.insert(GenomeId::new("genome2"));
        node.gene_members.insert(GenomeId::new("genome1"), vec!["geneA".to_string()]);
        node.gene_members.insert(GenomeId::new("genome2"), vec!["geneB".to_string()]);
        graph.add_node(node);

        GmlWriter::write(&graph, &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("length"), "GML should have length attribute");
        assert!(content.contains("seq"), "GML should have seq attribute");
        assert!(content.contains("protein"), "GML should have protein attribute");
        assert!(content.contains("genomes"), "GML should have genomes attribute");
        assert!(content.contains("member"), "GML should have member attribute");
    }
}
