//! P1 Integration tests for PanMiner.
//!
//! Tests P1 features:
//! - Contig-end pruning
//! - SV matrix output
//! - Paralog handling
//! - Large dataset support

use tempfile::TempDir;

use panminer::config::{PanminerConfig, CorrectionMode, OutputFormat};
use panminer::pipeline::PanminerPipeline;

// Test helpers for creating test data
mod test_helpers {
    use std::fs::File;
    use std::io::Write;
    use tempfile::TempDir;

    /// Create a minimal GFF3 file with gene annotations.
    ///
    /// The gene_id parameter allows specifying the same gene ID across multiple genomes
    /// so they will cluster together during pangenome analysis.
    pub fn create_test_gff(dir: &TempDir, name: &str, gene_id: &str, start: u32, end: u32) -> std::path::PathBuf {
        let gff_path = dir.path().join(format!("{}.gff", name));
        let mut file = File::create(&gff_path).unwrap();

        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "##sequence-region {} 1 10000", name).unwrap();
        writeln!(file, "{}\tProkka\tgene\t{}\t{}\t.\t+\t.\tID={};product=test_protein", name, start, end, gene_id).unwrap();

        writeln!(file, "##FASTA").unwrap();
        writeln!(file, ">{}", name).unwrap();
        // Write a sequence that's long enough to cover the gene
        let seq_len = (end + 100) as usize;
        let mut sequence = Vec::with_capacity(seq_len);
        for i in 0..seq_len {
            sequence.push(b"ATCG"[i % 4]);
        }
        let seq_str = String::from_utf8(sequence).unwrap();
        for i in (0..seq_str.len()).step_by(80) {
            writeln!(file, "{}", &seq_str[i..std::cmp::min(i + 80, seq_str.len())]).unwrap();
        }

        gff_path
    }

    /// Create a GFF with multiple genes on a single contig.
    /// This is used to test contig-end pruning with genes at contig boundaries.
    pub fn create_test_gff_with_contig_ends(dir: &TempDir, name: &str, gene_ids: &[&str]) -> std::path::PathBuf {
        let gff_path = dir.path().join(format!("{}.gff", name));
        let mut file = File::create(&gff_path).unwrap();

        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "##sequence-region {} 1 50000", name).unwrap();

        // Create genes at various positions on the contig
        // First gene (potential contig start)
        writeln!(file, "{}\tProkka\tgene\t100\t200\t.\t+\t.\tID={};product=test", name, gene_ids[0]).unwrap();
        // Middle genes
        for (i, &gene_id) in gene_ids.iter().enumerate().skip(1) {
            let start = 1000 + (i * 200) as u32;
            let end = start + 99;
            writeln!(file, "{}\tProkka\tgene\t{}\t{}\t.\t+\t.\tID={};product=test", name, start, end, gene_id).unwrap();
        }
        // Last gene (potential contig end)
        let end_gene_id = gene_ids[gene_ids.len() - 1];
        writeln!(file, "{}\tProkka\tgene\t49800\t49900\t.\t+\t.\tID={};product=test", name, end_gene_id).unwrap();

        writeln!(file, "##FASTA").unwrap();
        writeln!(file, ">{}", name).unwrap();
        // Write sequence (50000 bytes)
        let seq_len = 50000;
        let mut sequence = Vec::with_capacity(seq_len);
        for i in 0..seq_len {
            sequence.push(b"ATCG"[i % 4]);
        }
        let seq_str = String::from_utf8(sequence).unwrap();
        for i in (0..seq_str.len()).step_by(80) {
            writeln!(file, "{}", &seq_str[i..std::cmp::min(i + 80, seq_str.len())]).unwrap();
        }

        gff_path
    }

    /// Create a GFF with a paralog (duplicate gene in the same genome).
    pub fn create_test_gff_with_paralog(dir: &TempDir, name: &str) -> std::path::PathBuf {
        let gff_path = dir.path().join(format!("{}.gff", name));
        let mut file = File::create(&gff_path).unwrap();

        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "##sequence-region {} 1 20000", name).unwrap();

        // First copy of the gene
        writeln!(file, "{}\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene1;product=test", name).unwrap();
        // Second copy - a paralog (slightly different sequence, same name for testing)
        writeln!(file, "{}\tProkka\tgene\t10000\t10100\t.\t+\t.\tID=gene1;product=test_paralog", name).unwrap();

        writeln!(file, "##FASTA").unwrap();
        writeln!(file, ">{}", name).unwrap();
        // Write sequence (need to cover up to 10100)
        let seq_len = 10200;
        let mut sequence = Vec::with_capacity(seq_len);
        for i in 0..seq_len {
            sequence.push(b"ATCG"[i % 4]);
        }
        let seq_str = String::from_utf8(sequence).unwrap();
        for i in (0..seq_str.len()).step_by(80) {
            writeln!(file, "{}", &seq_str[i..std::cmp::min(i + 80, seq_str.len())]).unwrap();
        }

        gff_path
    }

    /// Create a GFF with a single gene on a contig (for testing contig-end cases).
    pub fn create_test_gff_single_gene(dir: &TempDir, name: &str, gene_id: &str) -> std::path::PathBuf {
        let gff_path = dir.path().join(format!("{}.gff", name));
        let mut file = File::create(&gff_path).unwrap();

        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "##sequence-region {} 1 5000", name).unwrap();
        writeln!(file, "{}\tProkka\tgene\t100\t200\t.\t+\t.\tID={};product=test", name, gene_id).unwrap();

        writeln!(file, "##FASTA").unwrap();
        writeln!(file, ">{}", name).unwrap();
        let seq_len = 300;
        let mut sequence = Vec::with_capacity(seq_len);
        for i in 0..seq_len {
            sequence.push(b"ATCG"[i % 4]);
        }
        let seq_str = String::from_utf8(sequence).unwrap();
        for i in (0..seq_str.len()).step_by(80) {
            writeln!(file, "{}", &seq_str[i..std::cmp::min(i + 80, seq_str.len())]).unwrap();
        }

        gff_path
    }
}

#[test]
fn test_pipeline_contig_end_pruning() {

    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create genomes with genes at contig ends
    // Genome1: has a single gene at position 100-200 (likely contig end)
    let gff1 = test_helpers::create_test_gff_single_gene(&temp_dir, "genome1", "gene1");
    // Genome2: has a single gene at a different position
    let gff2 = test_helpers::create_test_gff_single_gene(&temp_dir, "genome2", "gene1");

    // Configure pipeline with SV matrix output (uses same underlying graph structure)
    let config = PanminerConfig::new()
        .with_input_files(vec![gff1, gff2])
        .with_output_dir(output_dir.clone())
        .with_threads(2)
        .with_mode(CorrectionMode::Default)
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json, OutputFormat::SVMatrix].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);

    // Run pipeline - should complete without errors
    let result = pipeline.run();
    assert!(result.is_ok(), "Pipeline should complete successfully with contig-end genes");

    let output_paths = result.unwrap();

    // Verify output files were created
    assert!(output_paths.matrix_csv.is_some(), "CSV matrix should be created");
    assert!(output_paths.graph.is_some(), "GML graph should be created");
    assert!(output_paths.json.is_some(), "JSON summary should be created");
    assert!(output_paths.sv_matrix.is_some(), "SV matrix should be created");
}

#[test]
fn test_pipeline_sv_matrix_output() {
    use std::fs::read_to_string;

    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create multiple genomes with adjacent genes
    // This should produce structural variants (inversions, duplications)
    let gff1 = test_helpers::create_test_gff_with_contig_ends(&temp_dir, "genome1", &["g1a", "g1b", "g1c", "g1d"]);
    let gff2 = test_helpers::create_test_gff_with_contig_ends(&temp_dir, "genome2", &["g2a", "g2b", "g2c", "g2d"]);
    let gff3 = test_helpers::create_test_gff_with_contig_ends(&temp_dir, "genome3", &["g3a", "g3b", "g3c", "g3d"]);

    // Configure pipeline with SV matrix output
    let config = PanminerConfig::new()
        .with_input_files(vec![gff1, gff2, gff3])
        .with_output_dir(output_dir.clone())
        .with_threads(2)
        .with_mode(CorrectionMode::Default)
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json, OutputFormat::SVMatrix].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);

    let result = pipeline.run();
    assert!(result.is_ok(), "Pipeline should complete successfully with SV matrix output");

    let output_paths = result.unwrap();

    // Verify SV matrix file is created
    assert!(output_paths.sv_matrix.is_some(), "SV matrix file should be created");

    // Verify the SV matrix file exists and has content
    let sv_path = output_paths.sv_matrix.unwrap();
    assert!(sv_path.exists(), "SV matrix file should exist");
    let content = read_to_string(&sv_path).unwrap();
    // The file should have a header and at least some data
    assert!(content.contains("VariantID"), "SV matrix should have VariantID header");
}

#[test]
fn test_pipeline_with_paralogs() {

    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create genomes with paralogs (duplicate genes)
    let gff1 = test_helpers::create_test_gff_with_paralog(&temp_dir, "genome1");
    let gff2 = test_helpers::create_test_gff_with_paralog(&temp_dir, "genome2");

    // Configure pipeline
    let config = PanminerConfig::new()
        .with_input_files(vec![gff1, gff2])
        .with_output_dir(output_dir.clone())
        .with_threads(2)
        .with_mode(CorrectionMode::Default)
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);

    // Run pipeline - should handle paralogs without errors
    let result = pipeline.run();
    assert!(result.is_ok(), "Pipeline should handle paralogs");

    let output_paths = result.unwrap();

    // Verify outputs were created
    assert!(output_paths.matrix_csv.is_some(), "CSV matrix should be created");
    assert!(output_paths.json.is_some(), "JSON summary should be created");
}

#[test]
fn test_pipeline_large_dataset() {
    use std::fs::read_to_string;

    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create 10+ genomes for a large dataset test
    let mut gff_paths = Vec::new();
    for i in 1..=12 {
        let gff = test_helpers::create_test_gff(&temp_dir, &format!("large_genome{}", i), "common_gene", 100, 200);
        gff_paths.push(gff);
    }

    // Configure pipeline
    let config = PanminerConfig::new()
        .with_input_files(gff_paths)
        .with_output_dir(output_dir.clone())
        .with_threads(2)
        .with_mode(CorrectionMode::Default)
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);

    // Run pipeline - should handle large dataset without errors
    let result = pipeline.run();
    assert!(result.is_ok(), "Pipeline should handle large dataset with 10+ genomes");

    let output_paths = result.unwrap();

    // Verify all outputs were created
    assert!(output_paths.matrix_csv.is_some(), "CSV matrix should be created");
    assert!(output_paths.graph.is_some(), "GML graph should be created");
    assert!(output_paths.json.is_some(), "JSON summary should be created");

    // Verify matrix contains all 12 genomes
    let matrix_path = output_paths.matrix_csv.unwrap();
    let content = read_to_string(&matrix_path).unwrap();
    for i in 1..=12 {
        assert!(content.contains(&format!("large_genome{}", i)), "Matrix should contain large_genome{}", i);
    }
}
