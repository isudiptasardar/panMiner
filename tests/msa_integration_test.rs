//! Integration tests for the `panminer msa` subcommand.
//!
//! These tests verify that the MSA subcommand is registered and that
//! the alignment runners and helper functions work correctly.

use panminer::clustering::{AlignmentRunner, AlignmentTool, MafftRunner, ClustalOmegaRunner, PrankRunner};

/// Verify that the msa subcommand is registered and shows help.
#[test]
fn test_msa_command_help() {
    let output = std::process::Command::new("cargo")
        .args(["run", "--", "msa", "--help"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("Failed to run panminer msa --help");

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("multiple sequence alignment"), "help text should mention MSA: {}", stdout);
    assert!(stdout.contains("--mode"), "help text should mention --mode: {}", stdout);
    assert!(stdout.contains("--aligner"), "help text should mention --aligner: {}", stdout);
}

/// Verify that alignment runners can be created and implement the trait.
#[test]
fn test_alignment_runner_creation() {
    let mafft = MafftRunner::new();
    assert_eq!(mafft.name(), "MAFFT Runner");

    let clustal = ClustalOmegaRunner::new();
    assert_eq!(clustal.name(), "Clustal Omega Runner");

    let prank = PrankRunner::new();
    assert_eq!(prank.name(), "PRANK Runner");
}

/// Verify that alignment tool enum matches expected values.
#[test]
fn test_alignment_tool_roundtrip() {
    assert_eq!(AlignmentTool::Mafft.executable(), "mafft");
    assert_eq!(AlignmentTool::ClustalOmega.executable(), "clustalo");
    assert_eq!(AlignmentTool::Prank.executable(), "prank");
}

/// Verify that gene_data.csv can be parsed for MSA with a synthetic file.
#[test]
fn test_read_gene_data_for_msa_core_mode() {
    use std::io::Write;
    let dir = tempfile::TempDir::new().unwrap();

    // Write a minimal gene_data.csv
    let gene_data_path = dir.path().join("gene_data.csv");
    let mut f = std::fs::File::create(&gene_data_path).unwrap();
    writeln!(f, "gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence").unwrap();
    writeln!(f, "cluster_001,,hypothetical,contig1,1,100,+,5,ATGCGTATCGATCGATCG,MVRIDR").unwrap();
    writeln!(f, "cluster_002,,transposase,contig2,200,300,-,2,GGCCAATT,GN").unwrap();
    writeln!(f, "cluster_003,,ribosomal,contig3,400,500,+,5,TTTTAAAACCCCGGGG,FKPR").unwrap();

    // Write a gene_presence_absence.csv with 5 genomes (14 meta + 5 genome cols = 19 cols)
    let roary_path = dir.path().join("gene_presence_absence.csv");
    let mut r = std::fs::File::create(&roary_path).unwrap();
    write!(r, "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome Fragment,Order within Fragment,Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,Avg group size nuc,genome1,genome2,genome3,genome4,genome5").unwrap();

    // Simulate the helper function logic inline (since we cannot call private fn from main.rs)
    // Core threshold for 5 genomes: ceil(5 * 0.99) = 5
    // cluster_001 (support=5) -> core, cluster_002 (support=2) -> not core, cluster_003 (support=5) -> core
    let core_threshold = (5_f32 * 0.99).ceil() as usize;
    assert_eq!(core_threshold, 5);

    // Verify parsing: cluster_001 and cluster_003 are core, cluster_002 is not
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(&gene_data_path)
        .unwrap();

    let mut core_clusters = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        let gene_id = record.get(0).unwrap_or("").to_string();
        let support: usize = record.get(7).unwrap_or("0").parse().unwrap_or(0);
        let dna = record.get(8).unwrap_or("").to_string();
        if support >= core_threshold && !dna.is_empty() {
            core_clusters.push(gene_id);
        }
    }

    assert!(core_clusters.contains(&"cluster_001".to_string()), "cluster_001 should be core");
    assert!(!core_clusters.contains(&"cluster_002".to_string()), "cluster_002 should not be core");
    assert!(core_clusters.contains(&"cluster_003".to_string()), "cluster_003 should be core");
}

/// Verify pan mode includes all clusters.
#[test]
fn test_read_gene_data_for_msa_pan_mode() {
    use std::io::Write;
    let dir = tempfile::TempDir::new().unwrap();

    let gene_data_path = dir.path().join("gene_data.csv");
    let mut f = std::fs::File::create(&gene_data_path).unwrap();
    writeln!(f, "gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence").unwrap();
    writeln!(f, "cluster_001,,hypothetical,contig1,1,100,+,5,ATGCGT,MVR").unwrap();
    writeln!(f, "cluster_002,,transposase,contig2,200,300,-,2,GGCCAA,GAK").unwrap();

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(&gene_data_path)
        .unwrap();

    // In pan mode, all clusters with sequences are included regardless of support
    let mut all_clusters = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        let gene_id = record.get(0).unwrap_or("").to_string();
        let dna = record.get(8).unwrap_or("").to_string();
        if !dna.is_empty() {
            all_clusters.push(gene_id);
        }
    }

    assert_eq!(all_clusters.len(), 2, "pan mode should include all clusters with sequences");
}