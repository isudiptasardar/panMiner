//! PanMiner CLI entry point.

use std::path::PathBuf;
use clap::Parser;
use tracing_subscriber::EnvFilter;

use panminer::config::{CorrectionMode, OutputFormat, PanminerConfig, QcMode};
use panminer::pipeline::PanminerPipeline;

/// PanMiner - A modern pangenome analysis tool with GPU and CPU support.
#[derive(Parser, Debug)]
#[command(name = "panminer", version, about)]
struct Cli {
    /// Input GFF3 files
    #[arg(required = true)]
    input: Vec<PathBuf>,

    /// Output directory
    #[arg(short, long, default_value = "panminer_output")]
    output: PathBuf,

    /// Number of threads (0 = auto-detect)
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,

    /// Chunk size for streaming processing
    #[arg(long, default_value = "100")]
    chunk_size: usize,

    /// Clustering identity threshold (0.5-1.0)
    #[arg(long, default_value = "0.98")]
    identity: f32,

    /// Correction mode: strict, default, sensitive
    #[arg(long, default_value = "default")]
    mode: String,

    /// Force CPU processing (disable GPU)
    #[arg(long)]
    force_cpu: bool,

    /// Disable MMseqs2 clustering
    #[arg(long)]
    no_mmseqs2: bool,

    /// Disable GPU detection and acceleration
    #[arg(long)]
    no_gpu: bool,

    /// Path to MMseqs2 binary
    #[arg(long)]
    mmseqs_path: Option<PathBuf>,

    /// Disable pre-processing QC (Mash/CheckM)
    #[arg(long)]
    no_qc: bool,

    /// QC mode: strict, default, sensitive
    #[arg(long, default_value = "default")]
    qc_mode: String,

    /// Path to CheckM2 database
    #[arg(long)]
    checkm_database: Option<PathBuf>,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,

    /// Output formats (comma-separated: matrix,alignment,graph,json)
    #[arg(long, default_value = "matrix,alignment,graph")]
    formats: String,
}

fn parse_mode(s: &str) -> CorrectionMode {
    match s.to_lowercase().as_str() {
        "strict" => CorrectionMode::Strict,
        "sensitive" => CorrectionMode::Sensitive,
        _ => CorrectionMode::Default,
    }
}

fn parse_qc_mode(s: &str) -> QcMode {
    match s.to_lowercase().as_str() {
        "strict" => QcMode::Strict,
        "sensitive" => QcMode::Sensitive,
        _ => QcMode::Default,
    }
}

fn parse_formats(s: &str) -> std::collections::HashSet<OutputFormat> {
    s.split(',')
        .filter_map(|f| match f.trim().to_lowercase().as_str() {
            "matrix" => Some(OutputFormat::Matrix),
            "alignment" => Some(OutputFormat::Alignment),
            "graph" => Some(OutputFormat::Graph),
            "json" => Some(OutputFormat::Json),
            "parquet" => Some(OutputFormat::Parquet),
            "html" => Some(OutputFormat::HtmlViz),
            _ => None,
        })
        .collect()
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    // Initialize logging
    let filter = if cli.verbose {
        EnvFilter::new("debug")
    } else {
        EnvFilter::new("info")
    };

    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(false)
        .init();

    tracing::info!("PanMiner v{}", panminer::VERSION);

    // Build configuration
    let config = PanminerConfig::new()
        .with_input_files(cli.input)
        .with_output_dir(cli.output)
        .with_threads(cli.threads)
        .with_chunk_size(cli.chunk_size)
        .with_mode(parse_mode(&cli.mode))
        .with_outputs(parse_formats(&cli.formats))
        .force_cpu(cli.force_cpu)
        .with_enable_mmseqs(!cli.no_mmseqs2)
        .with_prefer_gpu(!cli.no_gpu)
        .with_enable_qc(!cli.no_qc)
        .with_qc_mode(parse_qc_mode(&cli.qc_mode));

    // Validate
    config.validate()?;

    tracing::info!(
        "Processing {} genomes with {} threads",
        config.input_files.len(),
        config.effective_threads()
    );

    // Run pipeline
    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run()?;

    tracing::info!("Output written to: {:?}", result.output_dir);
    tracing::info!("Done.");

    Ok(())
}
