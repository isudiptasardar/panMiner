//! Core/accessory alignment output (FASTA format).

use std::io::Write;
use std::path::Path;

use crate::error::Result;
use crate::graph::PangenomeGraph;

/// Writer for alignment output in FASTA format.
pub struct AlignmentWriter;

impl AlignmentWriter {
    /// Write core gene alignments to a FASTA file.
    ///
    /// Core genes are those present in all (or nearly all) genomes.
    pub fn write_core(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        let total_genomes = graph.genomes.len().max(1);

        // Find core clusters (present in >= 99% of genomes)
        let core_threshold = (total_genomes as f32 * 0.99).ceil() as usize;

        for (cluster_id, node) in &graph.nodes {
            if node.support >= core_threshold {
                // Write cluster as a FASTA entry
                let annotation = node.annotations
                    .iter()
                    .next()
                    .cloned()
                    .unwrap_or_else(|| "hypothetical protein".to_string());

                writeln!(file, ">{} {}", cluster_id, annotation)?;
                // Placeholder: in full implementation, write aligned sequences
                writeln!(file, "# Cluster support: {}/{}", node.support, total_genomes)?;
            }
        }

        Ok(())
    }

    /// Write accessory gene alignments to a FASTA file.
    ///
    /// Accessory genes are those present in some but not all genomes.
    pub fn write_accessory(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        let total_genomes = graph.genomes.len().max(1);
        let core_threshold = (total_genomes as f32 * 0.99).ceil() as usize;

        for (cluster_id, node) in &graph.nodes {
            if node.support > 0 && node.support < core_threshold {
                let annotation = node.annotations
                    .iter()
                    .next()
                    .cloned()
                    .unwrap_or_else(|| "hypothetical protein".to_string());

                writeln!(file, ">{} {} support={}/{}", cluster_id, annotation, node.support, total_genomes)?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_core_empty_graph() {
        let graph = PangenomeGraph::new();
        let temp = NamedTempFile::new().unwrap();
        AlignmentWriter::write_core(&graph, temp.path()).unwrap();

        let content = std::fs::read_to_string(temp.path()).unwrap();
        assert!(content.is_empty());
    }
}
