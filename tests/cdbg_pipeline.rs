//! Integration tests for the cDBG-based pipeline mode.
//!
//! These tests verify the cDBG pipeline infrastructure without requiring
//! GGCAT or ggCaller to be installed (they test config, error handling,
//! and module wiring).

use panminer::config::{PanminerConfig, PipelineMode};
use panminer::error::Error;

#[test]
fn test_pipeline_mode_dbg_config() {
    let config = PanminerConfig::new()
        .with_pipeline_mode(PipelineMode::Dbg)
        .with_kmer_size(31);
    assert_eq!(config.pipeline_mode, PipelineMode::Dbg);
    assert_eq!(config.kmer_size, 31);
}

#[test]
fn test_pipeline_mode_gff_config() {
    let config = PanminerConfig::new()
        .with_pipeline_mode(PipelineMode::Gff);
    assert_eq!(config.pipeline_mode, PipelineMode::Gff);
}

#[test]
fn test_pipeline_mode_default_is_gff() {
    let config = PanminerConfig::default();
    assert_eq!(config.pipeline_mode, PipelineMode::Gff);
}

#[test]
fn test_pipeline_mode_from_str() {
    use std::str::FromStr;
    assert_eq!(PipelineMode::from_str("dbg").unwrap(), PipelineMode::Dbg);
    assert_eq!(PipelineMode::from_str("gff").unwrap(), PipelineMode::Gff);
    assert_eq!(PipelineMode::from_str("DBG").unwrap(), PipelineMode::Dbg);
    assert!(PipelineMode::from_str("invalid").is_err());
}

#[test]
fn test_dbg_mode_kmer_validation() {
    // k-mer too small
    let config = PanminerConfig::new()
        .with_pipeline_mode(PipelineMode::Dbg)
        .with_kmer_size(10)
        .with_input_files(vec!["test.fa".into()]);
    let err = config.validate().unwrap_err();
    assert!(err.to_string().contains("k-mer size must be >= 15"), "Expected k-mer error, got: {}", err);

    // k-mer too large
    let config = PanminerConfig::new()
        .with_pipeline_mode(PipelineMode::Dbg)
        .with_kmer_size(200)
        .with_input_files(vec!["test.fa".into()]);
    let err = config.validate().unwrap_err();
    assert!(err.to_string().contains("k-mer size must be <= 127"), "Expected k-mer error, got: {}", err);
}

#[test]
fn test_ggcaller_runner_detect() {
    // Verify detect() doesn't panic (may or may not find ggCaller)
    let _ = panminer::io::GGCallerRunner::detect();
}

#[test]
fn test_ggcaller_runner_path() {
    let runner = panminer::io::GGCallerRunner::new(
        std::path::PathBuf::from("/usr/local/bin/ggcaller"),
    );
    assert_eq!(runner.path(), std::path::Path::new("/usr/local/bin/ggcaller"));
}

#[test]
fn test_ggcaller_not_found_error() {
    let err = Error::ggcaller_not_found();
    assert!(err.to_string().contains("ggCaller not found"));
}

#[test]
fn test_ggcat_build_failed_error() {
    let err = Error::ggcat_build_failed("out of memory");
    assert!(err.to_string().contains("GGCAT cDBG build failed"));
}

#[test]
fn test_dbg_mode_errors_without_feature() {
    // When the dbg feature is not enabled, running the pipeline with
    // PipelineMode::Dbg should return FeatureNotEnabled, not silently continue.
    #[cfg(not(feature = "dbg"))]
    {
        use panminer::pipeline::PanminerPipeline;
        use std::path::PathBuf;

        let config = PanminerConfig::new()
            .with_input_files(vec![PathBuf::from("nonexistent.gff")])
            .with_output_dir(PathBuf::from("test_cdbg_no_feature"))
            .with_pipeline_mode(PipelineMode::Dbg);

        let pipeline = PanminerPipeline::new(config);
        let result = pipeline.run();
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("dbg") || msg.contains("feature"),
            "Expected feature-related error, got: {}",
            msg
        );
    }
}