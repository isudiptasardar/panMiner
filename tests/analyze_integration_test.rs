//! Integration tests for the `panminer analyze` subcommand downstream tool wiring.
//!
//! These tests verify that all downstream runners can be created, detected,
//! and implement the DownstreamRunner trait correctly.

use std::path::PathBuf;
use panminer::downstream::{DownstreamRunner, DownstreamInput};

/// Verify that Scoary2Runner is detected (or not) without panicking.
#[test]
fn test_scoary2_detect_no_panic() {
    let _ = panminer::downstream::Scoary2Runner::detect();
}

/// Verify that SpydrPickRunner is detected (or not) without panicking.
#[test]
fn test_spydrpick_detect_no_panic() {
    let _ = panminer::downstream::SpydrPickRunner::detect();
}

/// Verify that PanstripeRunner is detected (or not) without panicking.
#[test]
fn test_panstripe_detect_no_panic() {
    let _ = panminer::downstream::PanstripeRunner::detect();
}

/// Verify that AmrFinderRunner is detected (or not) without panicking.
#[test]
fn test_amrfinder_detect_no_panic() {
    let _ = panminer::downstream::AmrFinderRunner::detect();
}

/// Verify that Scoary2Runner implements DownstreamRunner trait.
#[test]
fn test_scoary2_implements_downstream_runner() {
    if let Some(runner) = panminer::downstream::Scoary2Runner::detect() {
        assert_eq!(runner.name(), "Scoary2");
        assert!(runner.is_available());
        let inputs = runner.required_inputs();
        assert!(inputs.contains(&DownstreamInput::PresenceAbsenceCsv));
    }
}

/// Verify that SpydrPickRunner implements DownstreamRunner trait.
#[test]
fn test_spydrpick_implements_downstream_runner() {
    if let Some(runner) = panminer::downstream::SpydrPickRunner::detect() {
        assert_eq!(runner.name(), "SpydrPick");
        assert!(runner.is_available());
        let inputs = runner.required_inputs();
        assert!(inputs.contains(&DownstreamInput::PresenceAbsenceCsv));
    }
}

/// Verify that PanstripeRunner implements DownstreamRunner trait.
#[test]
fn test_panstripe_implements_downstream_runner() {
    if let Some(runner) = panminer::downstream::PanstripeRunner::detect() {
        assert_eq!(runner.name(), "Panstripe");
        assert!(runner.is_available());
        let inputs = runner.required_inputs();
        assert!(inputs.contains(&DownstreamInput::PresenceAbsenceCsv));
        assert!(inputs.contains(&DownstreamInput::PhylogeneticTree));
    }
}

/// Verify that AmrFinderRunner implements DownstreamRunner trait.
#[test]
fn test_amrfinder_implements_downstream_runner() {
    if let Some(runner) = panminer::downstream::AmrFinderRunner::detect() {
        assert_eq!(runner.name(), "AMRFinderPlus");
        assert!(runner.is_available());
        let inputs = runner.required_inputs();
        assert!(inputs.contains(&DownstreamInput::ProteinFasta));
        assert!(inputs.contains(&DownstreamInput::GeneDataCsv));
    }
}

/// Test that all downstream runners can be created with builder methods.
#[test]
fn test_downstream_runner_builder_methods() {
    let _ = panminer::downstream::Scoary2Runner::new()
        .with_phenotypes(PathBuf::from("phenos.txt"))
        .with_output_dir(PathBuf::from("output"))
        .with_threads(4);

    let _ = panminer::downstream::SpydrPickRunner::new()
        .with_output_dir(PathBuf::from("output"))
        .with_threads(4);

    let _ = panminer::downstream::PanstripeRunner::new()
        .with_tree(PathBuf::from("tree.nwk"))
        .with_output_dir(PathBuf::from("output"));

    let _ = panminer::downstream::AmrFinderRunner::new()
        .with_database(PathBuf::from("/path/to/db"))
        .with_organism("Escherichia coli".to_string())
        .with_threads(4);
}