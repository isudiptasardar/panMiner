//! Benchmark suite for PanMiner.
//!
//! This module provides performance benchmarks for the main pipeline phases:
//! - GFF3 parsing
//! - Graph construction
//! - Error correction
//! - Output generation

use criterion::{criterion_group, criterion_main, Criterion};
use tempfile::TempDir;

use panminer::config::{PanminerConfig, OutputFormat};
use panminer::pipeline::PanminerPipeline;

/// Create a test GFF3 file with the specified number of genomes and genes per genome.
fn create_test_gff(dir: &TempDir, genome_name: &str, gene_count: usize) -> std::path::PathBuf {
    use std::fs::File;
    use std::io::Write;

    let gff_path = dir.path().join(format!("{}.gff", genome_name));
    let mut file = File::create(&gff_path).unwrap();

    writeln!(file, "##gff-version 3").unwrap();
    writeln!(file, "##sequence-region {} 1 {}", genome_name, gene_count * 200).unwrap();

    // Create genes at regular intervals
    for i in 0..gene_count {
        let start = 100 + i * 200;
        let end = start + 99;
        writeln!(
            file,
            "{}\tProkka\tgene\t{}\t{}\t.\t+\t.\tID=gene{};product=test_protein",
            genome_name, start, end, i
        )
        .unwrap();
    }

    writeln!(file, "##FASTA").unwrap();
    writeln!(file, ">{}", genome_name).unwrap();

    // Write sequence (gene_count * 200 bytes, padded)
    let seq_len = gene_count * 200 + 100;
    let sequence = "ATCG".repeat(seq_len / 4 + 1);
    for i in (0..seq_len).step_by(80) {
        writeln!(
            file,
            "{}",
            &sequence[i..std::cmp::min(i + 80, seq_len)]
        )
        .unwrap();
    }

    gff_path
}

/// Benchmark full pipeline with small dataset.
fn benchmark_small_pipeline(c: &mut Criterion) {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create 5 genomes with 100 genes each
    let mut gff_paths = Vec::new();
    for i in 0..5 {
        let gff = create_test_gff(&temp_dir, &format!("genome{}", i), 100);
        gff_paths.push(gff);
    }

    c.bench_function("pipeline_5_genomes_100_genes", |b| {
        b.iter(|| {
            let config = PanminerConfig::new()
                .with_input_files(gff_paths.clone())
                .with_output_dir(output_dir.clone())
                .with_threads(2)
                .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

            let pipeline = PanminerPipeline::new(config);
            let _result = pipeline.run();
        })
    });
}

/// Benchmark full pipeline with medium dataset.
fn benchmark_medium_pipeline(c: &mut Criterion) {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create 10 genomes with 200 genes each
    let mut gff_paths = Vec::new();
    for i in 0..10 {
        let gff = create_test_gff(&temp_dir, &format!("medium_genome{}", i), 200);
        gff_paths.push(gff);
    }

    c.bench_function("pipeline_10_genomes_200_genes", |b| {
        b.iter(|| {
            let config = PanminerConfig::new()
                .with_input_files(gff_paths.clone())
                .with_output_dir(output_dir.clone())
                .with_threads(4)
                .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

            let pipeline = PanminerPipeline::new(config);
            let _result = pipeline.run();
        })
    });
}

/// Benchmark full pipeline with large dataset.
fn benchmark_large_pipeline(c: &mut Criterion) {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create 25 genomes with 300 genes each
    let mut gff_paths = Vec::new();
    for i in 0..25 {
        let gff = create_test_gff(&temp_dir, &format!("large_genome{}", i), 300);
        gff_paths.push(gff);
    }

    c.bench_function("pipeline_25_genomes_300_genes", |b| {
        b.iter(|| {
            let config = PanminerConfig::new()
                .with_input_files(gff_paths.clone())
                .with_output_dir(output_dir.clone())
                .with_threads(8)
                .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::Json].into_iter().collect());

            let pipeline = PanminerPipeline::new(config);
            let _result = pipeline.run();
        })
    });
}

criterion_group!(
    benches,
    benchmark_small_pipeline,
    benchmark_medium_pipeline,
    benchmark_large_pipeline,
);

criterion_main!(benches);
