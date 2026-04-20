//! PanMiner CLI entry point.

use std::path::PathBuf;
use clap::{Parser, Subcommand};
use tracing_subscriber::EnvFilter;

use panminer::config::{CorrectionMode, OutputFormat, PanminerConfig, PipelineMode, QcMode};
use panminer::io::BaktaDbType;
use panminer::output::{filter_presence_absence, parse_filter_types};
use panminer::pipeline::PanminerPipeline;
use panminer::io::QcRunner;
use panminer::downstream::{DownstreamRunner, DownstreamResult};
use panminer::clustering::{AlignmentRunner, AlignmentTool, MafftRunner, ClustalOmegaRunner, PrankRunner};

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

    /// Length difference cutoff for clustering (0.0-1.0).
    /// Gene pairs with relative length difference > (1 - len_dif_percent) are excluded.
    #[arg(long, default_value = "0.98")]
    len_dif_percent: f32,

    /// Correction mode: strict, default, sensitive
    #[arg(long, default_value = "default")]
    mode: String,

    /// Pipeline mode: gff (GFF3-annotated) or dbg (cDBG-based gene calling)
    #[arg(long, default_value = "gff")]
    pipeline_mode: String,

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

    /// Output formats (comma-separated: matrix,alignment,graph,gff,json)
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

    /// Alignment filtering method: none, clipkit, bmge
    #[arg(long, default_value = "none")]
    filter_alignment: String,

    /// Run GWAS analysis with pyseer after pangenome construction
    #[arg(long)]
    gwas: bool,

    /// Path to phenotype file for GWAS (TSV: genome_id phenotype_value)
    #[arg(long)]
    phenotype: Option<PathBuf>,

    /// k-mer size for cDBG construction (only used with --pipeline-mode dbg)
    #[arg(long, default_value = "31")]
    kmer_size: usize,

    /// Collapsing thresholds (comma-separated, high to low)
    #[arg(long, default_value = "0.99,0.95,0.9,0.8,0.7")]
    collapse_thresholds: String,
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

        /// Generate MDS scatter plot and HTML report
        #[arg(long)]
        mds: bool,

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

    /// Extract gene sequences for a specific cluster
    #[command(name = "extract-gene")]
    ExtractGene {
        /// PanMiner output directory
        #[arg(short = 'i', long)]
        input: PathBuf,

        /// Cluster ID to extract
        #[arg(long)]
        cluster: String,

        /// Output FASTA file
        #[arg(short = 'o', long, default_value = "extracted_genes.fasta")]
        output: PathBuf,

        /// Extract protein sequences instead of DNA
        #[arg(long)]
        protein: bool,
    },

    /// Integrate a new genome into an existing pangenome
    #[command(name = "integrate")]
    Integrate {
        /// Existing PanMiner output directory
        #[arg(long)]
        graph: PathBuf,

        /// New GFF3 file to integrate
        #[arg(long)]
        input: PathBuf,

        /// Output directory for updated pangenome
        #[arg(short = 'o', long, default_value = "integrated_output")]
        output: PathBuf,

        /// Identity threshold for matching new genes (0.5-1.0)
        #[arg(long, default_value = "0.98")]
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

        /// Generate gene abundance visualization (HTML/D3.js)
        #[arg(long)]
        abundance: bool,

        /// Run Pangrowth pangenome openness estimation
        #[arg(long)]
        pangrowth: bool,

        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },

    /// Run multiple sequence alignment on pangenome output
    #[command(name = "msa")]
    Msa {
        /// PanMiner output directory
        #[arg(short = 'i', long)]
        input: PathBuf,

        /// Output directory for alignment files
        #[arg(short = 'o', long, default_value = "msa_output")]
        output: PathBuf,

        /// Alignment mode: core or pan
        #[arg(long, default_value = "core")]
        mode: String,

        /// Alignment tool: mafft, clustal, prank
        #[arg(long, default_value = "mafft")]
        aligner: String,

        /// Number of threads
        #[arg(short = 't', long, default_value = "1")]
        threads: usize,

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
            "gff" => Some(OutputFormat::Gff),
            "json" => Some(OutputFormat::Json),
            "parquet" => Some(OutputFormat::Parquet),
            "html" => Some(OutputFormat::HtmlViz),
            "struct" => Some(OutputFormat::Struct),
            "svmatrix" => Some(OutputFormat::SVMatrix),
            "abundance" => Some(OutputFormat::Abundance),
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

fn parse_pipeline_mode(s: &str) -> PipelineMode {
    match s.to_lowercase().as_str() {
        "dbg" => PipelineMode::Dbg,
        "prodigal" => PipelineMode::Prodigal,
        _ => PipelineMode::Gff,
    }
}

fn parse_filter_method(s: &str) -> panminer::config::FilterMethod {
    match s.to_lowercase().as_str() {
        "bmge" => panminer::config::FilterMethod::Bmge,
        "clipkit" => panminer::config::FilterMethod::ClipKit,
        _ => panminer::config::FilterMethod::None,
    }
}

/// Count total genomes from the gene_presence_absence.csv header.
///
/// The Roary-compatible CSV has 14 metadata columns; everything after that
/// is a per-genome column. Returns 0 if the file is missing or unparseable.
fn count_genomes_from_roary_csv(input_dir: &PathBuf) -> usize {
    let roary_path = input_dir.join("gene_presence_absence.csv");
    if !roary_path.exists() {
        tracing::debug!("gene_presence_absence.csv not found, defaulting genome count to 0");
        return 0;
    }

    match std::fs::read_to_string(&roary_path) {
        Ok(content) => {
            let header = match content.lines().next() {
                Some(h) => h,
                None => return 0,
            };
            let num_columns = header.split(',').count();
            num_columns.saturating_sub(14)
        }
        Err(_) => 0,
    }
}

/// A row from gene_data.csv relevant for MSA.
struct GeneDataMsaRow {
    gene_id: String,
    support: usize,
    dna_sequence: String,
}

/// Read gene_data.csv and return clusters filtered by mode.
///
/// Returns a Vec of (cluster_id, sequences) where each cluster has one
/// centroid sequence. When `mode` is "core", only clusters present in
/// >= 99% of genomes are included. When "pan", all clusters are included.
fn read_gene_data_for_msa(
    path: &PathBuf,
    mode: &str,
    total_genomes: usize,
) -> anyhow::Result<Vec<(String, Vec<(String, Vec<u8>)>)>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)?;

    let headers = reader.headers()?.clone();
    let column_map = build_column_map(&headers);

    let mut rows: Vec<GeneDataMsaRow> = Vec::new();

    for result in reader.records() {
        let record = result?;
        let gene_id = get_field(&record, &column_map, "gene_id").to_string();
        let support: usize = get_field(&record, &column_map, "support")
            .parse()
            .unwrap_or(0);
        let dna_sequence = get_field(&record, &column_map, "dna_sequence").to_string();

        rows.push(GeneDataMsaRow {
            gene_id,
            support,
            dna_sequence,
        });
    }

    let core_threshold = if total_genomes > 0 {
        (total_genomes as f32 * 0.99).ceil() as usize
    } else {
        // Fallback: if we could not determine genome count from Roary CSV,
        // try to infer from the max support value. Use the max support as
        // total_genomes so that only clusters with support == max are "core".
        let max_support = rows.iter().map(|r| r.support).max().unwrap_or(0);
        (max_support as f32 * 0.99).ceil() as usize
    };

    let clusters: Vec<(String, Vec<(String, Vec<u8>)>)> = rows
        .into_iter()
        .filter(|row| {
            if mode == "core" {
                row.support >= core_threshold
            } else {
                true // "pan" mode includes all
            }
        })
        .filter_map(|row| {
            if row.dna_sequence.is_empty() {
                return None;
            }
            let gene_id = row.gene_id.clone();
            Some((
                gene_id,
                vec![(row.gene_id, row.dna_sequence.into_bytes())],
            ))
        })
        .collect();

    Ok(clusters)
}

/// Build a mapping from column name to column index.
fn build_column_map(headers: &csv::StringRecord) -> std::collections::HashMap<String, usize> {
    headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect()
}

/// Get a field from a CSV record by column name, falling back to empty string.
fn get_field<'a>(
    record: &'a csv::StringRecord,
    column_map: &std::collections::HashMap<String, usize>,
    name: &str,
) -> &'a str {
    column_map
        .get(name)
        .and_then(|&i| record.get(i))
        .unwrap_or("")
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
        Some(Commands::Qc { input, output, qc_mode: _, distance, mds, verbose }) => {
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

            // Optionally compute distance matrix and MDS projection
            let mds_projection: Option<panminer::io::MdsProjection> = if (distance || mds) && !input.is_empty() {
                if let Some(dist) = qc_runner.compute_distance_matrix(&input) {
                    let dist_path = output.join("distance_matrix.csv");
                    dist.write_csv(&dist_path)?;
                    tracing::info!("Wrote distance matrix to: {:?}", dist_path);

                    // Compute MDS with real genome labels
                    let labels: Vec<String> = input.iter()
                        .filter_map(|p| p.file_stem().and_then(|s| s.to_str()).map(|s| s.to_string()))
                        .collect();
                    if let Ok(mds_proj) = panminer::io::compute_mds_with_labels(&dist.distance_matrix, &labels) {
                        let mds_path = output.join("mds_coordinates.csv");
                        mds_proj.write_csv(&mds_path)?;
                        tracing::info!("Wrote MDS coordinates to: {:?}", mds_path);
                        Some(mds_proj)
                    } else {
                        None
                    }
                } else {
                    tracing::warn!("No distance tool available. Install skani: conda install -c bioconda skani");
                    None
                }
            } else {
                None
            };

            // Write HTML report if requested
            if mds {
                let html_path = output.join("qc_report.html");
                panminer::output::write_qc_html_report(&qc_results, mds_projection.as_ref(), &html_path)?;
                tracing::info!("Wrote QC HTML report to: {:?}", html_path);
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
        Some(Commands::ExtractGene { input, cluster, output, protein }) => {
            tracing::info!("PanMiner extract-gene v{}", panminer::VERSION);
            panminer::io::extract_gene(&input, &cluster, &output, protein)?;
            tracing::info!("Extracted sequences written to: {:?}", output);
            tracing::info!("Done.");
        }
        Some(Commands::Integrate { graph, input, output, identity, threads, verbose }) => {
            let filter = if verbose {
                EnvFilter::new("debug")
            } else {
                EnvFilter::new("info")
            };
            tracing_subscriber::fmt()
                .with_env_filter(filter)
                .with_target(false)
                .init();

            tracing::info!("PanMiner integrate v{}", panminer::VERSION);
            let effective_threads = if threads == 0 {
                std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1)
            } else {
                threads
            };

            match panminer::io::integrate_genome(&graph, &input, &output, identity, effective_threads) {
                Ok(result) => {
                    tracing::info!("Integrated genome into pangenome");
                    tracing::info!("Total nodes: {}, Total edges: {}", result.total_nodes, result.total_edges);
                    tracing::info!("Genes matched to existing clusters: {}", result.genes_matched);
                    tracing::info!("New clusters created: {}", result.new_nodes);
                    tracing::info!("Output written to: {:?}", result.output_dir);
                }
                Err(e) => {
                    tracing::error!("Integration failed: {}", e);
                    return Err(anyhow::anyhow!("{}", e));
                }
            }
            tracing::info!("Done.");
        }
        Some(Commands::Analyze { input, gwas, gwas_tool, phenotypes, panstripe, tree, amr, amr_database, organism, neighborhood, seed_gene, neighborhood_depth, accumulation, num_samples, export_grapetree, export_itol, abundance, pangrowth, verbose }) => {
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
                        if panminer::downstream::PyseerRunner::is_installed() {
                            if phenotypes.is_none() {
                                tracing::warn!("pyseer requires --phenotypes file. Provide with: --phenotypes phenotypes.txt");
                            } else {
                                let runner = panminer::downstream::PyseerRunner::new()
                                    .with_phenotypes(phenotypes.as_ref().unwrap().clone());
                                match runner.run(&input) {
                                    Ok(result) => {
                                        result.write_to(&downstream_dir)?;
                                        tracing::info!("Pyseer analysis complete: {}", result.summary());
                                    }
                                    Err(e) => {
                                        tracing::error!("Pyseer analysis failed: {}", e);
                                    }
                                }
                            }
                        } else {
                            tracing::warn!("pyseer is not installed. Install with: pip install pyseer");
                        }
                    }
                    "scoary2" => {
                        if let Some(runner) = panminer::downstream::Scoary2Runner::detect() {
                            let mut runner = runner.with_output_dir(downstream_dir.clone());
                            if let Some(ref phenotypes_file) = phenotypes {
                                runner = runner.with_phenotypes(phenotypes_file.clone());
                            }
                            match runner.run(&input) {
                                Ok(result) => {
                                    result.write_to(&downstream_dir)?;
                                    tracing::info!("Scoary2 analysis complete: {}", result.summary());
                                }
                                Err(e) => {
                                    tracing::error!("Scoary2 analysis failed: {}", e);
                                }
                            }
                        } else {
                            tracing::warn!("scoary2 is not installed. Install with: pip install scoary-2");
                        }
                    }
                    "spydrpick" => {
                        if let Some(runner) = panminer::downstream::SpydrPickRunner::detect() {
                            let runner = runner.with_output_dir(downstream_dir.clone());
                            match runner.run(&input) {
                                Ok(result) => {
                                    result.write_to(&downstream_dir)?;
                                    tracing::info!("SpydrPick analysis complete: {}", result.summary());
                                }
                                Err(e) => {
                                    tracing::error!("SpydrPick analysis failed: {}", e);
                                }
                            }
                        } else {
                            tracing::warn!("spydrpick is not installed. Install with: conda install -c bioconda spydrpick");
                        }
                    }
                    _ => {
                        tracing::warn!("Unknown GWAS tool: {}. Available: pyseer, scoary2, spydrpick", gwas_tool);
                    }
                }
            }

            // Panstripe evolutionary model
            if panstripe {
                if let Some(runner) = panminer::downstream::PanstripeRunner::detect() {
                    let mut runner = runner.with_output_dir(downstream_dir.clone());
                    if let Some(ref tree_file) = tree {
                        runner = runner.with_tree(tree_file.clone());
                    }
                    match runner.run(&input) {
                        Ok(result) => {
                            result.write_to(&downstream_dir)?;
                            tracing::info!("Panstripe analysis complete: {}", result.summary());
                        }
                        Err(e) => {
                            tracing::error!("Panstripe analysis failed: {}", e);
                        }
                    }
                } else {
                    tracing::warn!("panstripe is not installed. Install with: conda install -c conda-forge r-panstripe");
                }
            }

            // AMRFinderPlus resistome analysis
            if amr {
                if let Some(runner) = panminer::downstream::AmrFinderRunner::detect() {
                    let mut runner = runner;
                    if let Some(ref db_path) = amr_database {
                        runner = runner.with_database(db_path.clone());
                    }
                    if let Some(ref org) = organism {
                        runner = runner.with_organism(org.clone());
                    }
                    match runner.run(&input) {
                        Ok(result) => {
                            result.write_to(&downstream_dir)?;
                            tracing::info!("AMRFinderPlus analysis complete: {}", result.summary());
                        }
                        Err(e) => {
                            tracing::error!("AMRFinderPlus analysis failed: {}", e);
                        }
                    }
                } else {
                    tracing::warn!("amrfinder is not installed. Install with: conda install -c bioconda ncbi-amrfinder");
                }
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

            // Pangrowth pangenome openness
            if pangrowth {
                if let Some(runner) = panminer::downstream::PangrowthRunner::detect() {
                    let runner = runner.with_output_dir(downstream_dir.clone());
                    match runner.run(&input) {
                        Ok(result) => {
                            result.write_to(&downstream_dir)?;
                            tracing::info!("Pangrowth analysis complete: {}", result.summary());
                        }
                        Err(e) => {
                            tracing::error!("Pangrowth analysis failed: {}", e);
                        }
                    }
                } else {
                    tracing::warn!("pangrowth is not installed. Install with: conda install -c bioconda pangrowth");
                }
            }

            // Gene abundance visualization
            if abundance {
                tracing::info!("Generating gene abundance visualization");

                let rtab_path = input.join("gene_presence_absence.Rtab");
                if rtab_path.exists() {
                    match std::fs::read_to_string(&rtab_path) {
                        Ok(content) => {
                            let mut lines = content.lines();
                            let header = lines.next().unwrap_or("");
                            let genome_names: Vec<&str> = header
                                .split('\t')
                                .skip(1) // skip "Gene" column
                                .collect();
                            let total_genomes = genome_names.len();

                            // Compute presence counts: histogram of (n_genomes, n_gene_families)
                            let mut count_histogram: std::collections::HashMap<usize, usize> =
                                std::collections::HashMap::new();
                            let mut total_clusters: usize = 0;

                            // For rarefaction: track which genomes each gene appears in
                            let mut gene_genome_sets: Vec<Vec<usize>> = Vec::new();

                            for line in lines {
                                let fields: Vec<&str> = line.split('\t').collect();
                                if fields.len() < 2 {
                                    continue;
                                }
                                let mut present_in: Vec<usize> = Vec::new();
                                for (i, val) in fields.iter().skip(1).enumerate() {
                                    if let Ok(v) = val.parse::<u8>() {
                                        if v > 0 {
                                            present_in.push(i);
                                        }
                                    }
                                }
                                let n = present_in.len();
                                *count_histogram.entry(n).or_insert(0) += 1;
                                total_clusters += 1;
                                gene_genome_sets.push(present_in);
                            }

                            // Build sorted presence_counts
                            let mut presence_counts: Vec<(usize, usize)> =
                                count_histogram.into_iter().collect();
                            presence_counts.sort_by_key(|&(n, _)| n);

                            // Compute rarefaction data
                            // Use a deterministic ordering: iterate over genome indices
                            // and count cumulative unique genes as genomes are added
                            let genome_order: Vec<usize> = (0..total_genomes).collect();
                            let mut rarefaction_data: Vec<(usize, usize)> = Vec::new();
                            let mut seen_genes: std::collections::HashSet<usize> =
                                std::collections::HashSet::new();

                            for (step, &genome_idx) in genome_order.iter().enumerate() {
                                for (gene_idx, genomes) in gene_genome_sets.iter().enumerate() {
                                    if genomes.contains(&genome_idx) {
                                        seen_genes.insert(gene_idx);
                                    }
                                }
                                rarefaction_data.push((step + 1, seen_genes.len()));
                            }

                            let output_path = downstream_dir.join("abundance_report.html");
                            match panminer::output::AbundanceVizWriter::write_report(
                                &presence_counts,
                                &rarefaction_data,
                                total_genomes,
                                total_clusters,
                                &output_path,
                            ) {
                                Ok(_) => {
                                    tracing::info!("Abundance report written to: {:?}", output_path);
                                }
                                Err(e) => {
                                    tracing::error!("Failed to write abundance report: {}", e);
                                }
                            }
                        }
                        Err(e) => {
                            tracing::error!("Failed to read Rtab file: {}", e);
                        }
                    }
                } else {
                    tracing::warn!(
                        "gene_presence_absence.Rtab not found in {:?}. Run PanMiner with --formats matrix first.",
                        input
                    );
                }
            }

            tracing::info!("Downstream analysis complete. Results in: {:?}", downstream_dir);
            tracing::info!("Done.");
        }
        Some(Commands::Msa { input, output, mode, aligner, threads, verbose }) => {
            let filter = if verbose {
                EnvFilter::new("debug")
            } else {
                EnvFilter::new("info")
            };
            tracing_subscriber::fmt()
                .with_env_filter(filter)
                .with_target(false)
                .init();

            tracing::info!("PanMiner msa v{}", panminer::VERSION);

            // Validate input directory
            let gene_data_path = input.join("gene_data.csv");
            if !gene_data_path.exists() {
                return Err(anyhow::anyhow!(
                    "gene_data.csv not found in {:?}. Is this a PanMiner output directory?",
                    input
                ));
            }

            // Create output directory
            std::fs::create_dir_all(&output)?;

            // Determine total genome count from gene_presence_absence.csv header
            let total_genomes = count_genomes_from_roary_csv(&input);

            // Build the alignment runner based on --aligner
            let (runner, tool): (Box<dyn AlignmentRunner>, AlignmentTool) = match aligner.as_str() {
                "mafft" => (Box::new(MafftRunner::new()), AlignmentTool::Mafft),
                "clustal" => (Box::new(ClustalOmegaRunner::new()), AlignmentTool::ClustalOmega),
                "prank" => (Box::new(PrankRunner::new()), AlignmentTool::Prank),
                other => {
                    return Err(anyhow::anyhow!(
                        "Unknown aligner '{}'. Available: mafft, clustal, prank",
                        other
                    ));
                }
            };

            if !runner.is_available() {
                return Err(anyhow::anyhow!(
                    "{} is not installed. Install with: conda install -c bioconda {}",
                    tool.name(),
                    tool.executable()
                ));
            }

            // Read gene_data.csv and collect sequences based on mode
            let clusters = read_gene_data_for_msa(&gene_data_path, &mode, total_genomes)?;

            if clusters.is_empty() {
                tracing::warn!("No gene clusters match the '{}' filter", mode);
                tracing::info!("Done.");
                return Ok(());
            }

            tracing::info!(
                "Running MSA on {} clusters with {} (mode: {})",
                clusters.len(),
                tool.name(),
                mode
            );

            let effective_threads = if threads == 0 {
                std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1)
            } else {
                threads
            };

            let mut aligned = 0usize;
            let mut failed = 0usize;

            for (cluster_id, sequences) in &clusters {
                if sequences.len() < 2 {
                    tracing::debug!("Skipping cluster '{}': only {} sequence(s)", cluster_id, sequences.len());
                    continue;
                }

                match runner.run_msa(sequences, tool) {
                    Ok(result) => {
                        let out_path = output.join(format!("{}.aligned.fasta", cluster_id));
                        std::fs::write(&out_path, &result.aligned_fasta)?;
                        tracing::debug!(
                            "Aligned cluster '{}': {} sequences, {} columns",
                            cluster_id, result.num_sequences, result.alignment_length
                        );
                        aligned += 1;
                    }
                    Err(e) => {
                        tracing::warn!("Alignment failed for cluster '{}': {}", cluster_id, e);
                        failed += 1;
                    }
                }

                // Simple concurrency hint for the OS
                if aligned % effective_threads == 0 && aligned > 0 {
                    std::thread::yield_now();
                }
            }

            tracing::info!(
                "MSA complete: {} clusters aligned, {} failed, {} skipped (single sequence)",
                aligned,
                failed,
                clusters.len() - aligned - failed
            );
            tracing::info!("Output written to: {:?}", output);
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
                .with_len_dif_percent(cli.len_dif_percent)
                .with_enable_qc(!cli.no_qc)
                .with_qc_mode(parse_qc_mode(&cli.qc_mode))
                .with_reannotate(cli.reannotate)
                .with_keep_bakta_output(cli.keep_bakta_output)
                .with_no_bakta_db_download(cli.no_bakta_db_download)
                .with_trim_alignment(cli.trim_alignment)
                .with_trim_mode(cli.trim_mode)
                .with_codons(cli.codons)
                .with_run_gwas(cli.gwas)
                .with_filter_method(parse_filter_method(&cli.filter_alignment));

            // Parse collapsing thresholds (comma-separated, high to low)
            let collapse_thresholds: Vec<f32> = cli.collapse_thresholds
                .split(',')
                .filter_map(|s| s.trim().parse().ok())
                .collect();
            if !collapse_thresholds.is_empty() {
                config = config.with_collapse_thresholds(collapse_thresholds);
            }

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

            config = config
                .with_pipeline_mode(parse_pipeline_mode(&cli.pipeline_mode))
                .with_kmer_size(cli.kmer_size);

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