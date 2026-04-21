//! Core/accessory alignment output (FASTA format) with MSA integration.
//!
//! Aligns each gene cluster separately (one MSA per cluster), then
//! concatenates the per-gene alignments into core and accessory files.
//! This matches Panaroo's approach where each gene family is aligned
//! independently, producing meaningful alignments for phylogenetics.

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

    /// Compute the core threshold (99% of genomes, at least 1).
    fn core_threshold(total_genomes: usize) -> usize {
        if total_genomes == 0 {
            return 1;
        }
        ((total_genomes as f64 * 0.99).ceil() as usize).max(1)
    }

    /// Write core gene alignments to a FASTA file.
    ///
    /// Core genes are those present in >= 99% of genomes.
    /// Each gene cluster is aligned separately (per-gene MSA), then
    /// all aligned sequences are concatenated into one output file.
    pub fn write_core(&self, graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        let total_genomes = graph.genomes.len();
        let core_threshold = Self::core_threshold(total_genomes);

        // Collect core clusters
        let core_clusters: Vec<_> = graph
            .nodes
            .iter()
            .filter(|(_, node)| node.support >= core_threshold)
            .collect();

        if core_clusters.is_empty() {
            writeln!(file, "# No core genes found")?;
            return Ok(());
        }

        // Align each core cluster separately
        let mut total_aligned = 0usize;
        let mut total_failed = 0usize;

        for (cluster_id, node) in &core_clusters {
            // Collect per-genome sequences for this cluster
            let sequences = Self::cluster_sequences(node, graph);

            if sequences.len() < 2 {
                // Single sequence: write unaligned
                if let Some((name, seq)) = sequences.first() {
                    writeln!(file, ">{}", name)?;
                    let seq_str = String::from_utf8_lossy(seq);
                    for chunk in seq_str.as_bytes().chunks(80) {
                        writeln!(file, "{}", String::from_utf8_lossy(chunk))?;
                    }
                }
                continue;
            }

            match self.runner.run_msa(&sequences, AlignmentTool::Mafft) {
                Ok(result) => {
                    file.write_all(result.aligned_fasta.as_bytes())?;
                    total_aligned += 1;
                }
                Err(e) => {
                    tracing::warn!("MSA failed for cluster {}: {}. Writing unaligned.", cluster_id, e);
                    // Fallback: write unaligned sequences
                    for (name, seq) in &sequences {
                        writeln!(file, ">{}", name)?;
                        let seq_str = String::from_utf8_lossy(seq);
                        for chunk in seq_str.as_bytes().chunks(80) {
                            writeln!(file, "{}", String::from_utf8_lossy(chunk))?;
                        }
                    }
                    total_failed += 1;
                }
            }
        }

        writeln!(file, "# Core genes (>=99% presence), Total genomes: {}", total_genomes)?;
        writeln!(file, "# Aligned clusters: {}, Failed: {}", total_aligned, total_failed)?;

        Ok(())
    }

    /// Write accessory gene alignments to a FASTA file.
    ///
    /// Accessory genes are those present in some but not all genomes.
    /// Each gene cluster is aligned separately.
    pub fn write_accessory(&self, graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        let total_genomes = graph.genomes.len();
        let core_threshold = Self::core_threshold(total_genomes);

        // Collect accessory clusters
        let acc_clusters: Vec<_> = graph
            .nodes
            .iter()
            .filter(|(_, node)| node.support > 0 && node.support < core_threshold)
            .collect();

        if acc_clusters.is_empty() {
            writeln!(file, "# No accessory genes found")?;
            return Ok(());
        }

        let mut total_aligned = 0usize;
        let mut total_failed = 0usize;

        for (cluster_id, node) in &acc_clusters {
            let sequences = Self::cluster_sequences(node, graph);

            if sequences.len() < 2 {
                if let Some((name, seq)) = sequences.first() {
                    writeln!(file, ">{}", name)?;
                    let seq_str = String::from_utf8_lossy(seq);
                    for chunk in seq_str.as_bytes().chunks(80) {
                        writeln!(file, "{}", String::from_utf8_lossy(chunk))?;
                    }
                }
                continue;
            }

            match self.runner.run_msa(&sequences, AlignmentTool::Mafft) {
                Ok(result) => {
                    file.write_all(result.aligned_fasta.as_bytes())?;
                    total_aligned += 1;
                }
                Err(e) => {
                    tracing::warn!("MSA failed for accessory cluster {}: {}. Writing unaligned.", cluster_id, e);
                    for (name, seq) in &sequences {
                        writeln!(file, ">{}", name)?;
                        let seq_str = String::from_utf8_lossy(seq);
                        for chunk in seq_str.as_bytes().chunks(80) {
                            writeln!(file, "{}", String::from_utf8_lossy(chunk))?;
                        }
                    }
                    total_failed += 1;
                }
            }
        }

        writeln!(file, "# Accessory genes (<99% presence), Total genomes: {}", total_genomes)?;
        writeln!(file, "# Aligned clusters: {}, Failed: {}", total_aligned, total_failed)?;

        Ok(())
    }

    /// Collect per-genome sequences for a cluster node.
    ///
    /// Each genome that has this cluster contributes one sequence.
    /// Falls back to centroid sequence for all genomes if per-genome
    /// sequences are not available.
    fn cluster_sequences(node: &crate::graph::Node, graph: &PangenomeGraph) -> Vec<(String, Vec<u8>)> {
        let mut sequences = Vec::new();
        let mut seen_seqs: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();

        for genome_id in &node.genomes {
            // Try to get per-genome gene sequence from gene_lookup
            let seq = if let Some(gene_ids) = node.gene_members.get(genome_id) {
                if let Some(gene_id_str) = gene_ids.first() {
                    if let Some(gene) = graph.gene_lookup.get(&crate::graph::GeneId::new(gene_id_str.as_str())) {
                        if !gene.sequence.is_empty() {
                            gene.sequence.clone()
                        } else {
                            continue;
                        }
                    } else {
                        continue;
                    }
                } else {
                    continue;
                }
            } else {
                // Fallback: use centroid sequence
                node.centroid_sequences.first().cloned().unwrap_or_default()
            };

            if seen_seqs.contains(&seq) {
                // Skip duplicate sequences — MSA on identical input is pointless
                continue;
            }
            seen_seqs.insert(seq.clone());

            sequences.push((format!("{}__{}", node.cluster_id, genome_id), seq));
        }

        // If no per-genome sequences found, use centroid once
        if sequences.is_empty() {
            if let Some(centroid) = node.centroid_sequences.first() {
                sequences.push((node.cluster_id.to_string(), centroid.clone()));
            }
        }

        sequences
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
    fn test_core_threshold() {
        assert_eq!(AlignmentWriter::core_threshold(0), 1);
        assert_eq!(AlignmentWriter::core_threshold(100), 99);
        assert_eq!(AlignmentWriter::core_threshold(1), 1);
        assert_eq!(AlignmentWriter::core_threshold(200), 198);
    }
}