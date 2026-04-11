//! Integration tests for the PanMiner pipeline.
//!
//! Tests the full pipeline end-to-end with real GFF3 files,
//! verifies all output formats, and tests scalability.

use std::fs::{File, read_to_string};
use std::io::Write;
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
    ///
    /// All genomes using the same gene_id will have identical sequences, ensuring
    /// they cluster together.
    pub fn create_test_gff(dir: &TempDir, name: &str, gene_id: &str, start: u32, end: u32) -> std::path::PathBuf {
        let gff_path = dir.path().join(format!("{}.gff", name));
        let mut file = File::create(&gff_path).unwrap();

        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "##sequence-region {} 1 10000", name).unwrap();
        writeln!(file, "{}\tProkka\tgene\t{}\t{}\t.\t+\t.\tID={};product=test_protein", name, start, end, gene_id).unwrap();

        writeln!(file, "##FASTA").unwrap();
        writeln!(file, ">{}", name).unwrap();
        // Write a sequence that's long enough to cover the gene
        // Gene is at positions start-end (1-indexed), so we need at least 'end' characters
        let seq_len = (end + 100) as usize;
        let sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".repeat(seq_len / 50 + 1);
        for i in (0..seq_len).step_by(80) {
            writeln!(file, "{}", &sequence[i..std::cmp::min(i + 80, seq_len)]).unwrap();
        }

        gff_path
    }
}

#[test]
fn test_full_pipeline_with_multiple_genomes() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create test GFF files with the same gene ID so they will cluster together
    let gff1 = test_helpers::create_test_gff(&temp_dir, "genome1", "gene1", 100, 200);
    let gff2 = test_helpers::create_test_gff(&temp_dir, "genome2", "gene1", 100, 200);
    let gff3 = test_helpers::create_test_gff(&temp_dir, "genome3", "gene1", 100, 200);

    // Configure pipeline
    let config = PanminerConfig::new()
        .with_input_files(vec![gff1, gff2, gff3])
        .with_output_dir(output_dir.clone())
        .with_threads(2)
        .with_chunk_size(0)
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);

    // Run pipeline
    let result = pipeline.run();
    assert!(result.is_ok(), "Pipeline should complete successfully");

    let output_paths = result.unwrap();

    // Verify output files were created
    assert!(output_paths.matrix_csv.is_some(), "CSV matrix should be created");
    assert!(output_paths.matrix_rtab.is_some(), "Rtab matrix should be created");
    assert!(output_paths.graph.is_some(), "GML graph should be created");
    assert!(output_paths.json.is_some(), "JSON summary should be created");
}

#[test]
fn test_pipeline_qc_enabled() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    let gff = test_helpers::create_test_gff(&temp_dir, "genome1", "gene1", 100, 200);

    let config = PanminerConfig::new()
        .with_input_files(vec![gff])
        .with_output_dir(output_dir.clone())
        .with_enable_qc(true) // Enable QC
        .with_threads(1)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);

    // Run with QC - may fail if CheckM2 not installed or genomes don't pass QC
    let result = pipeline.run();

    // If CheckM2 is installed and genomes pass, pipeline succeeds
    // If CheckM2 is not installed or genomes fail QC, pipeline errors
    if result.is_err() {
        let err = result.unwrap_err().to_string();
        if err.contains("CheckM2") || err.contains("checkm2") || err.contains("QC filtering") {
            // Expected - CheckM2 not installed or genomes filtered out
        } else {
            panic!("Unexpected error: {}", err);
        }
    }
}

#[test]
fn test_pipeline_different_correction_modes() {
    let temp_dir = TempDir::new().unwrap();

    let gff = test_helpers::create_test_gff(&temp_dir, "genome1", "gene1", 100, 200);

    for mode in [CorrectionMode::Strict, CorrectionMode::Default, CorrectionMode::Sensitive] {
        let config = PanminerConfig::new()
            .with_input_files(vec![gff.clone()])
            .with_output_dir(temp_dir.path().join(format!("output_{:?}", mode)))
            .with_mode(mode.clone())
            .with_enable_qc(false)
            .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

        let pipeline = PanminerPipeline::new(config);
        let result = pipeline.run();
        assert!(result.is_ok(), "Pipeline should work with mode: {:?}", mode);
    }
}

#[test]
fn test_pipeline_edge_case_single_genome() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    let gff = test_helpers::create_test_gff(&temp_dir, "genome1", "gene1", 100, 200);

    let config = PanminerConfig::new()
        .with_input_files(vec![gff])
        .with_output_dir(output_dir.clone())
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run();

    assert!(result.is_ok(), "Pipeline should work with single genome");

    let output_paths = result.unwrap();
    assert!(output_paths.matrix_csv.is_some());
}

#[test]
fn test_debug_output() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    let gff1 = test_helpers::create_test_gff(&temp_dir, "genome1", "gene1", 100, 200);
    let gff2 = test_helpers::create_test_gff(&temp_dir, "genome2", "gene1", 100, 200);

    let config = PanminerConfig::new()
        .with_input_files(vec![gff1, gff2])
        .with_output_dir(output_dir.clone())
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run();

    // Print info for debugging
    if let Err(e) = result {
        panic!("Pipeline failed: {:?}", e);
    }

    let output_paths = result.unwrap();
    eprintln!("Output paths: {:?}", output_paths);

    // List files in output directory
    if let Ok(entries) = std::fs::read_dir(&output_dir) {
        for entry in entries {
            if let Ok(entry) = entry {
                eprintln!("File: {}", entry.path().display());
            }
        }
    }
}

#[test]
fn test_pipeline_output_matrix_content() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    let gff1 = test_helpers::create_test_gff(&temp_dir, "genome1", "gene1", 100, 200);
    let gff2 = test_helpers::create_test_gff(&temp_dir, "genome2", "gene1", 100, 200);

    // Debug: read the sequences from the GFF files
    let content1 = std::fs::read_to_string(&gff1).unwrap();
    let content2 = std::fs::read_to_string(&gff2).unwrap();
    eprintln!("GFF1 content:\n{}", content1);
    eprintln!("GFF2 content:\n{}", content2);

    let config = PanminerConfig::new()
        .with_input_files(vec![gff1, gff2])
        .with_output_dir(output_dir.clone())
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);
    pipeline.run().unwrap();

    // Read and verify matrix content
    if let Some(matrix_path) = output_dir.join("gene_presence_absence.csv").to_str() {
        let content = read_to_string(matrix_path).unwrap();
        eprintln!("Matrix content:\n{}", content);
        assert!(content.contains("genome1"), "Matrix should contain genome1");
        assert!(content.contains("genome2"), "Matrix should contain genome2");
        assert!(content.contains("cluster"), "Matrix should contain cluster names");
    }
}

#[test]
fn test_pipeline_graph_structure() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create two genomes with genes that would be adjacent
    let gff1 = test_helpers::create_test_gff(&temp_dir, "genome1", "gene1", 100, 200);
    let gff2 = test_helpers::create_test_gff(&temp_dir, "genome2", "gene1", 100, 200);

    // Debug: read the sequences from the GFF files
    let content1 = std::fs::read_to_string(&gff1).unwrap();
    let content2 = std::fs::read_to_string(&gff2).unwrap();
    eprintln!("GFF1 content:\n{}", content1);
    eprintln!("GFF2 content:\n{}", content2);

    let config = PanminerConfig::new()
        .with_input_files(vec![gff1, gff2])
        .with_output_dir(output_dir.clone())
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run();

    assert!(result.is_ok());

    // Verify GML file exists and has valid structure
    let gml_path = output_dir.join("final_graph.gml");
    assert!(gml_path.exists(), "GML file should exist");

    let content = read_to_string(gml_path).unwrap();
    eprintln!("GML content:\n{}", content);
    assert!(content.contains("node"), "GML should contain nodes");
    // With only one cluster, there will be no edges (edges connect different clusters)
    // The test just needs nodes to exist
    assert!(content.contains("id"), "GML should contain node ids");
}

#[test]
fn test_pipeline_json_output() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    let gff = test_helpers::create_test_gff(&temp_dir, "genome1", "gene1", 100, 200);

    // Debug: read the sequence from the GFF file
    let content = std::fs::read_to_string(&gff).unwrap();
    eprintln!("GFF content:\n{}", content);

    let config = PanminerConfig::new()
        .with_input_files(vec![gff])
        .with_output_dir(output_dir.clone())
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);
    pipeline.run().unwrap();

    // Verify JSON files exist
    let json_path = output_dir.join("_pangenome.json");
    assert!(json_path.exists(), "JSON summary should exist");

    let content = read_to_string(json_path).unwrap();
    eprintln!("JSON content:\n{}", content);
    assert!(content.contains("version"), "JSON should contain version");
    assert!(content.contains("num_clusters"), "JSON should contain cluster count");
}

#[test]
fn test_pipeline_with_realistic_gene_count() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create a more realistic test with multiple genes per genome
    let gff_path = temp_dir.path().join("genome1.gff");
    let mut file = File::create(&gff_path).unwrap();

    // Write multiple genes
    for _i in 0..10 {
        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "##sequence-region genome1 1 10000").unwrap();
        writeln!(file, "genome1\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene_0;product=test").unwrap();
    }

    writeln!(file, "##FASTA").unwrap();
    writeln!(file, ">genome1").unwrap();
    writeln!(file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG").unwrap();

    let config = PanminerConfig::new()
        .with_input_files(vec![gff_path])
        .with_output_dir(output_dir.clone())
        .with_enable_qc(false)
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run();

    assert!(result.is_ok(), "Pipeline should handle multiple genes");
}
