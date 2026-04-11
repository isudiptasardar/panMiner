//! Integration tests for Bakta re-annotation.
//!
//! These tests verify that the Bakta re-annotation pipeline works correctly.
//! Tests that require Bakta to be installed are behind the `bakta` feature flag.

use panminer::io::{BaktaRunner, BaktaDbType, is_gff_file, is_genbank_file, genbank_to_fasta};
use std::path::PathBuf;

#[test]
fn test_bakta_db_type_display() {
    assert_eq!(BaktaDbType::Full.to_string(), "full");
    assert_eq!(BaktaDbType::Light.to_string(), "light");
}

#[test]
fn test_bakta_detect() {
    // This test will pass whether Bakta is installed or not
    // It just verifies detect() doesn't panic
    let result = BaktaRunner::detect();
    if result.is_some() {
        println!("Bakta detected: {:?}", result.unwrap().name());
    } else {
        println!("Bakta not installed (expected in CI)");
    }
}

#[test]
fn test_gff_file_detection() {
    // Verify that GFF files are correctly identified
    let gff_path = PathBuf::from("genome.gff");
    assert!(is_gff_file(&gff_path));

    let gff3_path = PathBuf::from("genome.gff3");
    assert!(is_gff_file(&gff3_path));

    let fasta_path = PathBuf::from("genome.fasta");
    assert!(!is_gff_file(&fasta_path));

    let gbk_path = PathBuf::from("genome.gbk");
    assert!(!is_gff_file(&gbk_path));
}

#[test]
fn test_genbank_file_detection() {
    assert!(is_genbank_file(&PathBuf::from("genome.gb")));
    assert!(is_genbank_file(&PathBuf::from("genome.gbk")));
    assert!(is_genbank_file(&PathBuf::from("genome.gbff")));
    assert!(is_genbank_file(&PathBuf::from("genome.genbank")));
    assert!(!is_genbank_file(&PathBuf::from("genome.fasta")));
    assert!(!is_genbank_file(&PathBuf::from("genome.gff")));
}

#[test]
fn test_genbank_conversion() {
    let temp_dir = tempfile::tempdir().unwrap();
    let gbk_path = temp_dir.path().join("test.gb");

    let gbk_content = r#"LOCUS       test                48 bp    DNA     linear   BCT 01-JAN-2024
DEFINITION  Test genome for conversion.
ACCESSION   test
VERSION     test.1
ORIGIN
        1 atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg atcgatcgat cg
//
"#;
    std::fs::write(&gbk_path, gbk_content).unwrap();

    let fasta_path = genbank_to_fasta(&gbk_path).unwrap();
    let fasta_content = std::fs::read_to_string(&fasta_path).unwrap();

    assert!(fasta_content.starts_with(">test\n"));
    // Verify nucleotide content (stripped of line numbers and spaces)
    let seq = fasta_content.strip_prefix(">test\n").unwrap();
    assert!(!seq.is_empty());
    assert!(seq.chars().all(|c| "atcgATCG".contains(c)));
}

#[test]
fn test_genbank_conversion_empty_origin() {
    let temp_dir = tempfile::tempdir().unwrap();
    let gbk_path = temp_dir.path().join("empty.gb");
    std::fs::write(&gbk_path, "LOCUS test\nDEFINITION empty\n//").unwrap();

    let result = genbank_to_fasta(&gbk_path);
    assert!(result.is_err());
}

#[test]
fn test_bakta_runner_builder() {
    let runner = BaktaRunner::new(
        PathBuf::from("bakta"),
        PathBuf::from("/path/to/db"),
    )
    .with_threads(8)
    .with_keep_contig_headers(false);

    assert_eq!(runner.name(), "Bakta");
    assert_eq!(runner.name_path(), PathBuf::from("bakta"));
}

#[test]
fn test_bakta_config_integration() {
    use panminer::config::PanminerConfig;

    let config = PanminerConfig::new()
        .with_reannotate(true)
        .with_bakta_threads(4)
        .with_no_bakta_db_download(true)
        .with_keep_bakta_output(true);

    assert!(config.reannotate);
    assert_eq!(config.bakta_threads, 4);
    assert!(config.no_bakta_db_download);
    assert!(config.keep_bakta_output);
}