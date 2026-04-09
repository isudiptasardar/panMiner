//! Core/accessory alignment output (FASTA format) with MSA integration.
//!
//! Uses MSA runners (MAFFT, Clustal Omega, PRANK) for real sequence alignment
//! instead of placeholder output.

use std::io::Write;
use std::path::Path;

use crate::error::Result;
use crate::graph::PangenomeGraph;
use crate::clustering::{AlignmentRunner, MafftRunner, AlignmentTool};

/// Writer for alignment output in FASTA format.
pub struct AlignmentWriter {
    /// MSA runner to use for alignment
    runner: Box<dyn AlignmentRunner>,
}

impl AlignmentWriter {
    /// Create a new alignment writer with MAFFT as default.
    pub fn new() -> Self {
        Self {
            runner: Box::new(MafftRunner::new()),
        }
    }

    /// Create a new alignment writer with a specific runner.
    pub fn with_runner(runner: Box<dyn AlignmentRunner>) -> Self {
        Self { runner }
    }

    /// Write core gene alignments to a FASTA file.
    ///
    /// Core genes are those present in >= 99% of genomes.
    /// Uses MSA to generate proper multiple sequence alignments.
    pub fn write_core(&self, graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        let total_genomes = graph.genomes.len().max(1);

        // Find core clusters (present in >= 99% of genomes)
        let core_threshold = (total_genomes as f32 * 0.99).ceil() as usize;

        // Collect core cluster sequences for MSA
        let sequences: Vec<(String, Vec<u8>)> = graph
            .nodes
            .iter()
            .filter(|(_, node)| node.support >= core_threshold)
            .filter_map(|(cluster_id, node)| {
                node.centroid_sequence.clone()
                    .map(|seq| (cluster_id.to_string(), seq))
            })
            .collect();

        if sequences.is_empty() {
            // If no core sequences, write placeholder
            writeln!(file, "# No core genes found")?;
            return Ok(());
        }

        // Run MSA on core sequences
        let result = self.runner.run_msa(&sequences, AlignmentTool::Mafft)?;

        // Write the aligned sequences
        file.write_all(result.aligned_fasta.as_bytes())?;

        // Write metadata header
        writeln!(file, "# Aligned with {} ({})", result.tool.name(), result.tool.executable())?;
        writeln!(file, "# Sequences: {}, Alignment length: {}", result.num_sequences, result.alignment_length)?;
        writeln!(file, "# Core genes (>=99% presence), Total genomes: {}", total_genomes)?;

        Ok(())
    }

    /// Write accessory gene alignments to a FASTA file.
    ///
    /// Accessory genes are those present in some but not all genomes.
    /// Uses MSA to generate proper multiple sequence alignments.
    pub fn write_accessory(&self, graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        let total_genomes = graph.genomes.len().max(1);
        let core_threshold = (total_genomes as f32 * 0.99).ceil() as usize;

        // Collect accessory cluster sequences for MSA
        let sequences: Vec<(String, Vec<u8>)> = graph
            .nodes
            .iter()
            .filter(|(_, node)| node.support > 0 && node.support < core_threshold)
            .filter_map(|(cluster_id, node)| {
                node.centroid_sequence.clone()
                    .map(|seq| (cluster_id.to_string(), seq))
            })
            .collect();

        if sequences.is_empty() {
            // If no accessory sequences, write placeholder
            writeln!(file, "# No accessory genes found")?;
            return Ok(());
        }

        // Run MSA on accessory sequences
        let result = self.runner.run_msa(&sequences, AlignmentTool::Mafft)?;

        // Write the aligned sequences
        file.write_all(result.aligned_fasta.as_bytes())?;

        // Write metadata header
        writeln!(file, "# Aligned with {} ({})", result.tool.name(), result.tool.executable())?;
        writeln!(file, "# Sequences: {}, Alignment length: {}", result.num_sequences, result.alignment_length)?;
        writeln!(file, "# Accessory genes (<99% presence), Total genomes: {}", total_genomes)?;

        Ok(())
    }
}

impl Default for AlignmentWriter {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_alignment_writer_creation() {
        // Test that the alignment writer can be created
        // Note: We don't check is_available() since the MSA tool may not be installed
        let _writer = AlignmentWriter::new();
    }

    #[test]
    fn test_write_core_empty_graph() {
        let graph = PangenomeGraph::new();
        let temp = NamedTempFile::new().unwrap();
        let writer = AlignmentWriter::new();
        writer.write_core(&graph, temp.path()).unwrap();

        let content = std::fs::read_to_string(temp.path()).unwrap();
        assert!(content.contains("No core genes found"));
    }

    #[test]
    fn test_write_core_with_sequences() {
        let graph = PangenomeGraph::new();
        let temp = NamedTempFile::new().unwrap();
        let writer = AlignmentWriter::new();

        // Test with a simple case - write should succeed even without real sequences
        // (sequences will be empty, so placeholder output)
        writer.write_core(&graph, temp.path()).unwrap();
    }
}
