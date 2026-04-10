//! Structural variant matrix output in TSV format.
//!
//! Outputs a presence/absence matrix for structural variants.
//! Format: First row is genome names, subsequent rows are variants with 1/0 presence.

use std::path::Path;
use std::fs::File;
use std::io::Write;

use crate::error::Result;
use crate::graph::{PangenomeGraph, StructuralVariantDetector};

/// Writer for structural variant matrix output.
pub struct SVMatrixWriter {
    genome_names: Vec<String>,
}

impl SVMatrixWriter {
    /// Create a new SV matrix writer.
    pub fn new() -> Self {
        Self {
            genome_names: Vec::new(),
        }
    }

    /// Set genome names for the matrix.
    pub fn with_genomes(mut self, names: Vec<String>) -> Self {
        self.genome_names = names;
        self
    }

    /// Extract SV triplets from the graph.
    /// Returns Vec<(variant_id, variant_type, presence_vector)>
    pub fn extract_triplets(graph: &PangenomeGraph) -> Vec<(String, String, Vec<bool>)> {
        let detector = StructuralVariantDetector::new();
        let variants = detector.detect(graph);

        // Get all genome names from the graph
        let genome_names: Vec<String> = graph
            .genomes
            .keys()
            .map(|g| g.as_str().to_string())
            .collect();

        // Build a map of genome indices
        let genome_index: std::collections::HashMap<String, usize> = genome_names
            .iter()
            .enumerate()
            .map(|(i, name)| (name.clone(), i))
            .collect();

        let mut triplets: Vec<(String, String, Vec<bool>)> = Vec::new();

        for variant in &variants {
            let variant_id = match variant.variant_type {
                crate::graph::VariantType::Inversion => format!("INV_{}", triplets.len() + 1),
                crate::graph::VariantType::Duplication => format!("DUP_{}", triplets.len() + 1),
                crate::graph::VariantType::Translocation => format!("TRA_{}", triplets.len() + 1),
                crate::graph::VariantType::Deletion => format!("DEL_{}", triplets.len() + 1),
            };

            let variant_type = match variant.variant_type {
                crate::graph::VariantType::Inversion => "INVERSION",
                crate::graph::VariantType::Duplication => "DUPLICATION",
                crate::graph::VariantType::Translocation => "TRANSLOCATION",
                crate::graph::VariantType::Deletion => "DELETION",
            }
            .to_string();

            // Build presence vector for all genomes
            let mut presence: Vec<bool> = vec![false; genome_names.len()];
            for genome in &variant.affected_genomes {
                if let Some(&idx) = genome_index.get(genome) {
                    presence[idx] = true;
                }
            }

            triplets.push((variant_id, variant_type, presence));
        }

        triplets
    }

    /// Write the SV matrix to TSV file.
    pub fn write_tsv(&self, triplets: &[(String, String, Vec<bool>)], path: &Path) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = std::io::BufWriter::new(file);

        // Write header row with genome names
        let header = format!("VariantID,VariantType,{}", self.genome_names.join("\t"));
        writeln!(writer, "{}", header)?;

        // Write each variant as a row
        for (variant_id, variant_type, presence) in triplets {
            let presence_str: Vec<&str> = presence
                .iter()
                .map(|&p| if p { "1" } else { "0" })
                .collect();

            let row = format!(
                "{},{},{}",
                variant_id,
                variant_type,
                presence_str.join("\t")
            );
            writeln!(writer, "{}", row)?;
        }

        Ok(())
    }
}

impl Default for SVMatrixWriter {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{PangenomeGraph, GenomeId, ClusterId, Node, Edge};

    #[test]
    fn test_new_writer() {
        let writer = SVMatrixWriter::new();
        assert!(writer.genome_names.is_empty());
    }

    #[test]
    fn test_with_genomes() {
        let writer = SVMatrixWriter::new()
            .with_genomes(vec!["genome1".to_string(), "genome2".to_string()]);
        assert_eq!(writer.genome_names.len(), 2);
        assert_eq!(writer.genome_names[0], "genome1");
    }

    #[test]
    fn test_extract_triplets_empty_graph() {
        let graph = PangenomeGraph::new();
        let triplets = SVMatrixWriter::extract_triplets(&graph);
        assert!(triplets.is_empty());
    }

    #[test]
    fn test_extract_triplets_with_variants() {
        // Create a graph with a duplication (high degree node)
        let mut graph = PangenomeGraph::new();

        // Create a cluster with 3+ connections (duplication)
        let cluster1 = ClusterId::new("C1");
        let mut node1 = Node::from_cluster(&{
            let mut c = crate::graph::GeneCluster::new("C1");
            c.support = 5;
            c
        });
        node1.genomes.insert(GenomeId::new("genome1"));
        node1.genomes.insert(GenomeId::new("genome2"));
        graph.add_node(node1);

        // Create a second node
        let cluster2 = ClusterId::new("C2");
        let mut node2 = Node::from_cluster(&{
            let mut c = crate::graph::GeneCluster::new("C2");
            c.support = 3;
            c
        });
        node2.genomes.insert(GenomeId::new("genome1"));
        graph.add_node(node2);

        // Create a third node
        let cluster3 = ClusterId::new("C3");
        let mut node3 = Node::from_cluster(&{
            let mut c = crate::graph::GeneCluster::new("C3");
            c.support = 3;
            c
        });
        node3.genomes.insert(GenomeId::new("genome1"));
        graph.add_node(node3);

        // Create a fourth node
        let cluster4 = ClusterId::new("C4");
        let mut node4 = Node::from_cluster(&{
            let mut c = crate::graph::GeneCluster::new("C4");
            c.support = 3;
            c
        });
        node4.genomes.insert(GenomeId::new("genome1"));
        graph.add_node(node4);

        // Create edges to give C1 high degree (duplication signal)
        let mut edge1 = Edge::new(cluster1.clone(), cluster2.clone());
        edge1.add_genome(GenomeId::new("genome1"));
        graph.add_edge(edge1);

        let mut edge2 = Edge::new(cluster1.clone(), cluster3.clone());
        edge2.add_genome(GenomeId::new("genome1"));
        graph.add_edge(edge2);

        let mut edge3 = Edge::new(cluster1.clone(), cluster4.clone());
        edge3.add_genome(GenomeId::new("genome1"));
        graph.add_edge(edge3);

        let triplets = SVMatrixWriter::extract_triplets(&graph);
        // Should detect at least one duplication due to high degree
        assert!(triplets.len() >= 1);
        assert_eq!(triplets[0].0, "DUP_1");
        assert_eq!(triplets[0].1, "DUPLICATION");
        // All genomes should be present in the triplet (even if empty, it's a Vec<bool>)
    }

    #[test]
    fn test_write_tsv() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("sv_matrix.tsv");

        let triplets = vec![
            ("DUP_1".to_string(), "DUPLICATION".to_string(), vec![true, false]),
            ("INV_1".to_string(), "INVERSION".to_string(), vec![true, true]),
        ];

        let writer = SVMatrixWriter::new()
            .with_genomes(vec!["genome1".to_string(), "genome2".to_string()]);
        writer.write_tsv(&triplets, &path).unwrap();

        let content = std::fs::read_to_string(path).unwrap();
        // Header uses CSV for first two columns, then TSV for genomes
        assert!(content.contains("VariantID,VariantType,genome1\tgenome2"));
        assert!(content.contains("DUP_1,DUPLICATION,1\t0"));
        assert!(content.contains("INV_1,INVERSION,1\t1"));
    }
}
