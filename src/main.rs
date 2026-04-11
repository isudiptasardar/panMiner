//! PanMiner CLI entry point.

use std::path::PathBuf;
use clap::{Parser, Subcommand};
use tracing_subscriber::EnvFilter;

use panminer::config::{CorrectionMode, OutputFormat, PanminerConfig, QcMode};
use panminer::io::BaktaDbType;
use panminer::output::{filter_presence_absence, parse_filter_types};
use panminer::pipeline::PanminerPipeline;
use panminer::io::QcRunner;
use panminer::downstream::{DownstreamRunner, DownstreamResult};

/// PanMiner - A modern pangenome analysis tool with GPU and CPU support.
#[derive(Parser, Debug)]
#[command(name = "panminer", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    /// Input GFF3 files (for default run command)
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

    /// Re-annotate input genomes with Bakta before analysis
    #[arg(short = 'r', long)]
    reannotate: bool,

    /// Path to Bakta database directory
    #[arg(long)]
    bakta_db: Option<PathBuf>,

    /// Bakta database type for auto-download (full or light)
    #[arg(long, default_value = "full")]
    bakta_db_type: String,

    /// Number of threads for Bakta (default: same as pipeline)
    #[arg(long)]
    bakta_threads: Option<usize>,

    /// Fail if Bakta database not found instead of auto-downloading
    #[arg(long)]
    no_bakta_db_download: bool,

    /// Keep Bakta output files after pipeline completes
    #[arg(long)]
    keep_bakta_output: bool,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,

    /// Output formats (comma-separated: matrix,alignment,graph,json)
    #[arg(long, default_value = "matrix,alignment,graph")]
    formats: String,

    /// Trim alignment with ClipKIT after MSA
    #[arg(long)]
    trim_alignment: bool,

    /// ClipKIT trim mode: smart-gap, gappyout, strict
    #[arg(long, default_value = "smart-gap")]
    trim_mode: String,

    /// Generate codon alignments via MACSE
    #[arg(long)]
    codons: bool,

    /// Run GWAS analysis with pyseer after pangenome construction
    #[arg(long)]
    gwas: bool,

    /// Path to phenotype file for GWAS (TSV: genome_id phenotype_value)
    #[arg(long)]
    phenotype: Option<PathBuf>,
}

/// Subcommands for PanMiner.
#[derive(Subcommand, Debug)]
enum Commands {
    /// Filter a presence/absence matrix
    #[command(name = "filter-pa")]
    FilterPa {
        /// Input gene_presence_absence.csv file
        #[arg(long)]
        input: PathBuf,

        /// Output filtered file
        #[arg(long)]
        output: PathBuf,

        /// Filter types (comma-separated: frag,pseudo,length)
        #[arg(long, default_value = "frag,pseudo")]
        filter_type: String,

        /// Length outlier threshold (proportion deviation from mode)
        #[arg(long, default_value = "0.5")]
        length_threshold: f32,
    },

    /// Run genome QC (CheckM2 completeness/contamination + distance estimation)
    #[command(name = "qc")]
    Qc {
        /// Input genome files (GFF3/FASTA/GenBank)
        #[arg(required = true)]
        input: Vec<PathBuf>,

        /// Output directory for QC results
        #[arg(short, long, default_value = "qc_output")]
        output: PathBuf,

        /// QC mode: strict, default, sensitive
        #[arg(long, default_value = "default")]
        qc_mode: String,

        /// Compute pairwise ANI/distance matrix
        #[arg(long)]
        distance: bool,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },

    /// Merge multiple PanMiner output directories
    #[command(name = "merge")]
    Merge {
        /// Input directories to merge (at least 2)
        #[arg(required = true)]
        directories: Vec<PathBuf>,

        /// Output directory for merged pangenome
        #[arg(short, long, default_value = "merged_output")]
        output: PathBuf,

        /// Identity threshold for centroid clustering (0.5-1.0)
        #[arg(long, default_value = "0.95")]
        identity: f32,

        /// Number of threads (0 = auto-detect)
        #[arg(short = 't', long, default_value = "0")]
        threads: usize,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },

    /// Run downstream analyses on PanMiner output
    #[command(name = "analyze")]
    Analyze {
        /// PanMiner output directory from a prior run
        #[arg(short = 'i', long)]
        input: PathBuf,

        /// Run GWAS analysis
        #[arg(long)]
        gwas: bool,

        /// GWAS tool: pyseer (default), scoary2, spydrpick
        #[arg(long, default_value = "pyseer")]
        gwas_tool: String,

        /// Path to phenotypes file (TSV: genome_id<tab>phenotype)
        #[arg(long)]
        phenotypes: Option<PathBuf>,

        /// Run Panstripe evolutionary model
        #[arg(long)]
        panstripe: bool,

        /// Phylogenetic tree in Newick format
        #[arg(long)]
        tree: Option<PathBuf>,

        /// Run AMRFinderPlus resistome analysis
        #[arg(long)]
        amr: bool,

        /// AMRFinderPlus database path
        #[arg(long)]
        amr_database: Option<PathBuf>,

        /// Organism for taxon-specific AMR detection (e.g., "Escherichia coli")
        #[arg(long)]
        organism: Option<String>,

        /// Extract gene neighborhood
        #[arg(long)]
        neighborhood: bool,

        /// Seed gene/cluster ID for neighborhood extraction
        #[arg(long)]
        seed_gene: Option<String>,

        /// Maximum neighborhood depth (default: 5)
        #[arg(long)]
        neighborhood_depth: Option<usize>,

        /// Generate gene accumulation curves
        #[arg(long)]
        accumulation: bool,

        /// Number of rarefaction samples (default: 100)
        #[arg(long)]
        num_samples: Option<usize>,

        /// Export for GrapeTree
        #[arg(long)]
        export_grapetree: bool,

        /// Export for iTOL
        #[arg(long)]
        export_itol: bool,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },
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
            "struct" => Some(OutputFormat::Struct),
            "svmatrix" => Some(OutputFormat::SVMatrix),
            _ => None,
        })
        .collect()
}

fn parse_bakta_db_type(s: &str) -> BaktaDbType {
    match s.to_lowercase().as_str() {
        "light" => BaktaDbType::Light,
        _ => BaktaDbType::Full,
    }
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

    match cli.command {
        Some(Commands::FilterPa { input, output, filter_type, length_threshold }) => {
            tracing::info!("PanMiner filter-pa v{}", panminer::VERSION);
            let types = parse_filter_types(&filter_type);
            filter_presence_absence(&input, &output, &types, length_threshold)?;
            tracing::info!("Filtered output written to: {:?}", output);
            tracing::info!("Done.");
        }
        Some(Commands::Qc { input, output, qc_mode: _, distance, verbose }) => {
            let filter = if verbose {
                EnvFilter::new("debug")
            } else {
                EnvFilter::new("info")
            };
            tracing_subscriber::fmt()
                .with_env_filter(filter)
                .with_target(false)
                .init();

            tracing::info!("PanMiner qc v{}", panminer::VERSION);

            // Create output directory
            std::fs::create_dir_all(&output)?;

            // Run CheckM2 QC
            let qc_runner = panminer::io::CheckmQcRunner::new();
            let mut qc_results = Vec::new();

            for genome_path in &input {
                match qc_runner.run_qc(genome_path) {
                    Ok(qc) => {
                        tracing::info!("QC for {}: completeness={:.1}%, contamination={:.1}%",
                            qc.genome_id, qc.completeness, qc.contamination);
                        qc_results.push(qc);
                    }
                    Err(e) => {
                        tracing::warn!("QC failed for {:?}: {}", genome_path, e);
                    }
                }
            }

            // Write QC stats
            let qc_stats_path = output.join("qc_stats.csv");
            panminer::output::write_qc_stats(&qc_results, &qc_stats_path)?;
            tracing::info!("Wrote QC stats to: {:?}", qc_stats_path);

            let qc_summary_path = output.join("qc_summary.txt");
            panminer::output::write_qc_summary(&qc_results, &qc_summary_path)?;
            tracing::info!("Wrote QC summary to: {:?}", qc_summary_path);

            // Optionally compute distance matrix
            if distance && !input.is_empty() {
                if let Some(dist) = qc_runner.compute_distance_matrix(&input) {
                    let dist_path = output.join("distance_matrix.csv");
                    dist.write_csv(&dist_path)?;
                    tracing::info!("Wrote distance matrix to: {:?}", dist_path);
                } else {
                    tracing::warn!("No distance tool available. Install FastANI or enable sourmash feature.");
                }
            }

            tracing::info!("QC complete. Results in: {:?}", output);
            tracing::info!("Done.");
        }
        Some(Commands::Merge { directories, output, identity, threads, verbose }) => {
            let filter = if verbose {
                EnvFilter::new("debug")
            } else {
                EnvFilter::new("info")
            };
            tracing_subscriber::fmt()
                .with_env_filter(filter)
                .with_target(false)
                .init();

            tracing::info!("PanMiner merge v{}", panminer::VERSION);
            let effective_threads = if threads == 0 {
                std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1)
            } else {
                threads
            };

            match panminer::merge_pangenomes(&directories, &output, identity, effective_threads) {
                Ok(result) => {
                    tracing::info!("Merged {} input directories", result.num_inputs);
                    tracing::info!("Total nodes: {}, Total edges: {}", result.total_nodes, result.total_edges);
                    tracing::info!("Merged nodes: {}", result.merged_nodes);
                    tracing::info!("Output written to: {:?}", result.output_dir);
                }
                Err(e) => {
                    tracing::error!("Merge failed: {}", e);
                    return Err(anyhow::anyhow!("{}", e));
                }
            }
            tracing::info!("Done.");
        }
        Some(Commands::Analyze { input, gwas, gwas_tool, phenotypes, panstripe, tree, amr, amr_database, organism, neighborhood, seed_gene, neighborhood_depth, accumulation, num_samples, export_grapetree, export_itol, verbose }) => {
            let filter = if verbose {
                EnvFilter::new("debug")
            } else {
                EnvFilter::new("info")
            };
            tracing_subscriber::fmt()
                .with_env_filter(filter)
                .with_target(false)
                .init();

            tracing::info!("PanMiner analyze v{}", panminer::VERSION);

            // Create downstream output directory
            let downstream_dir = input.join("downstream");
            std::fs::create_dir_all(&downstream_dir)?;

            // GWAS analysis
            if gwas {
                tracing::info!("Running GWAS analysis with {}", gwas_tool);
                match gwas_tool.as_str() {
                    "pyseer" => {
                        if panminer::gwas::PyseerRunner::is_installed() {
                            // PyseerRunner works on in-memory graph/matrix, so we need to load first
                            // For now, just verify the input directory has required files
                            let presence_absence = input.join("gene_presence_absence.csv");
                            if presence_absence.exists() {
                                let mut runner = panminer::gwas::PyseerRunner::new();
                                if let Some(ref phenotypes_file) = phenotypes {
                                    runner.with_phenotypes(phenotypes_file.clone());
                                }
                                tracing::info!("PyseerRunner is available and ready");
                                tracing::info!("GWAS results would be written to: {:?}", downstream_dir.join("pyseer_associations.csv"));
                            } else {
                                tracing::warn!("gene_presence_absence.csv not found in input directory");
                            }
                        } else {
                            tracing::warn!("pyseer is not installed. Install with: pip install pyseer");
                        }
                    }
                    "scoary2" => {
                        tracing::info!("Scoary2Runner not yet implemented - coming in Phase 2");
                    }
                    "spydrpick" => {
                        tracing::info!("SpydrPickRunner not yet implemented - coming in Phase 2");
                    }
                    _ => {
                        tracing::warn!("Unknown GWAS tool: {}. Available: pyseer, scoary2, spydrpick", gwas_tool);
                    }
                }
            }

            // Panstripe evolutionary model
            if panstripe {
                tracing::info!("PanstripeRunner not yet implemented - coming in Phase 2");
            }

            // AMRFinderPlus resistome analysis
            if amr {
                tracing::info!("AmrFinderRunner not yet implemented - coming in Phase 2");
            }

            // Gene neighborhood extraction
            if neighborhood {
                if let Some(ref seed) = seed_gene {
                    let depth = neighborhood_depth.unwrap_or(5);
                    tracing::info!("Extracting gene neighborhood for seed '{}' with depth {}", seed, depth);

                    use panminer::downstream::exploration::neighborhood::GeneNeighborhoodExtractor;
                    use panminer::graph::ClusterId;

                    let extractor = GeneNeighborhoodExtractor::new(ClusterId::new(seed.clone()), depth);
                    match extractor.run(&input) {
                        Ok(result) => {
                            result.write_to(&downstream_dir)?;
                            tracing::info!("Neighborhood extraction complete");
                        }
                        Err(e) => {
                            tracing::error!("Neighborhood extraction failed: {}", e);
                        }
                    }
                } else {
                    tracing::warn!("--neighborhood requires --seed-gene");
                }
            }

            // Gene accumulation curves
            if accumulation {
                tracing::info!("Generating gene accumulation curves");

                use panminer::downstream::exploration::accumulation::AccumulationCurveRunner;

                let mut runner = AccumulationCurveRunner::new();
                if let Some(n) = num_samples {
                    runner = runner.with_num_samples(n);
                }

                match runner.run(&input) {
                    Ok(result) => {
                        result.write_to(&downstream_dir)?;
                        tracing::info!("Accumulation curve complete");
                    }
                    Err(e) => {
                        tracing::error!("Accumulation curve failed: {}", e);
                    }
                }
            }

            // GrapeTree export
            if export_grapetree {
                tracing::info!("Exporting for GrapeTree");

                use panminer::downstream::exploration::grapetree::GrapeTreeExportRunner;

                let runner = GrapeTreeExportRunner::new(export_itol);
                match runner.run(&input) {
                    Ok(result) => {
                        result.write_to(&downstream_dir)?;
                        tracing::info!("GrapeTree export complete");
                    }
                    Err(e) => {
                        tracing::error!("GrapeTree export failed: {}", e);
                    }
                }
            }

            tracing::info!("Downstream analysis complete. Results in: {:?}", downstream_dir);
            tracing::info!("Done.");
        }
        None => {
            // Default: run the main pipeline
            tracing::info!("PanMiner v{}", panminer::VERSION);

            let mut config = PanminerConfig::new()
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
                .with_qc_mode(parse_qc_mode(&cli.qc_mode))
                .with_reannotate(cli.reannotate)
                .with_keep_bakta_output(cli.keep_bakta_output)
                .with_no_bakta_db_download(cli.no_bakta_db_download)
                .with_trim_alignment(cli.trim_alignment)
                .with_trim_mode(cli.trim_mode)
                .with_codons(cli.codons)
                .with_run_gwas(cli.gwas);

            if let Some(phenotype) = cli.phenotype {
                config = config.with_phenotype_file(phenotype);
            }

            if let Some(bakta_db) = cli.bakta_db {
                config = config.with_bakta_db_path(bakta_db);
            }
            if let Some(bakta_threads) = cli.bakta_threads {
                config = config.with_bakta_threads(bakta_threads);
            }
            config = config.with_bakta_db_type(parse_bakta_db_type(&cli.bakta_db_type));

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
        }
    }

    Ok(())
}