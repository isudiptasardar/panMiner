//! Output module for PanMiner.
//!
//! Generates multiple output formats from the pangenome graph,
//! following Panaroo/Roary naming conventions for compatibility:
//! - Presence/absence matrix (CSV/TSV/Rtab)
//! - Core/accessory alignments (ALN/FASTA)
//! - GML graph format
//! - GFF3 per-genome corrected annotations
//! - Panaroo-style reference files (pan_genome_reference.fa, gene_data.csv)
//! - JSON/JSONL
//! - Parquet (optional)
//! - Structural variant matrix (optional)
//! - Interactive HTML visualization (optional)

mod matrix;
mod alignment;
mod graph;
mod gff;
mod json;
mod struct_csv;
mod sv_matrix;
mod summary;
mod parquet;
mod html_viz;
mod filter_pa;
mod trim;
mod codon;
pub mod qc_stats;
pub mod qc_viz;

pub use matrix::MatrixWriter;
pub use alignment::AlignmentWriter;
pub use graph::GmlWriter;
pub use gff::write_gff_files;
pub use json::JsonWriter;
pub use struct_csv::write_structural_variants;
pub use sv_matrix::SVMatrixWriter;
pub use summary::write_summary_stats;
pub use codon::MacseRunner;
pub use filter_pa::{FilterType, filter_presence_absence, parse_filter_types};
pub use trim::{ClipKitRunner, TrimMode, BmgeRunner};
pub use qc_stats::{write_qc_stats, write_qc_summary};
pub use qc_viz::write_qc_html_report;

#[cfg(feature = "parquet")]
pub use parquet::ParquetWriter;

#[cfg(feature = "viz")]
pub use html_viz::HtmlVizWriter;

use std::path::PathBuf;
use std::collections::HashMap;

use crate::config::{OutputFormat, PanminerConfig, FilterMethod};
use crate::error::Result;
use crate::graph::PangenomeGraph;
use crate::graph::BitPackedMatrix;

/// Writes all requested output formats in parallel.
pub struct OutputWriter {
    output_dir: PathBuf,
    formats: Vec<OutputFormat>,
    trim_alignment: bool,
    trim_mode: String,
    codons: bool,
    filter_method: FilterMethod,
}

impl OutputWriter {
    /// Create a new output writer from config.
    pub fn new(config: &PanminerConfig) -> Self {
        Self {
            output_dir: config.output_dir.clone(),
            formats: config.outputs.iter().cloned().collect(),
            trim_alignment: config.trim_alignment,
            trim_mode: config.trim_mode.clone(),
            codons: config.codons,
            filter_method: config.filter_method,
        }
    }

    /// Create output directory if needed.
    pub fn ensure_output_dir(&self) -> Result<()> {
        std::fs::create_dir_all(&self.output_dir)?;
        Ok(())
    }

    /// Write all requested output formats in parallel.
    pub fn write_all(
        &self,
        graph: &PangenomeGraph,
        matrix: &BitPackedMatrix,
        gene_members: &HashMap<String, HashMap<String, Vec<String>>>,
    ) -> Result<OutputPaths> {
        self.ensure_output_dir()?;

        let mut paths = OutputPaths {
            output_dir: self.output_dir.clone(),
            matrix_csv: None,
            matrix_rtab: None,
            matrix_roary_csv: None,
            alignment: None,
            codon_alignment: None,
            bmge_alignment: None,
            graph: None,
            gff: None,
            reference_fasta: None,
            gene_data: None,
            dna_fasta: None,
            protein_fasta: None,
            json: None,
            jsonl: None,
            struct_csv: None,
            sv_matrix: None,
            summary_stats: None,
            pre_filt_graph: None,
            #[cfg(feature = "parquet")]
            parquet: None,
            #[cfg(feature = "viz")]
            html_viz: None,
        };

        // Write formats (sequential for now to collect paths)
        for format in &self.formats {
            match format {
                OutputFormat::Matrix => {
                    // Panaroo/Roary compatible CSV (14 metadata columns)
                    let csv_path = self.output_dir.join("gene_presence_absence.csv");
                    MatrixWriter::write_roary_csv(matrix, &csv_path)?;
                    paths.matrix_csv = Some(csv_path);
                    tracing::info!("Wrote gene presence/absence matrix (Roary-compatible CSV)");

                    // Roary-compatible Rtab (binary TSV)
                    let rtab_path = self.output_dir.join("gene_presence_absence.Rtab");
                    MatrixWriter::write_tsv(matrix, &rtab_path)?;
                    paths.matrix_rtab = Some(rtab_path);
                    tracing::info!("Wrote gene presence/absence matrix (Rtab)");

                    // Roary-compatible gene member CSV
                    let roary_gene_path = self.output_dir.join("gene_presence_absence_roary.csv");
                    match MatrixWriter::write_roary_gene_csv(matrix, gene_members, &roary_gene_path) {
                        Ok(_) => {
                            paths.matrix_roary_csv = Some(roary_gene_path);
                            tracing::info!("Wrote gene presence/absence Roary CSV with gene member IDs");
                        }
                        Err(e) => tracing::warn!("Failed to write Roary gene CSV: {}", e),
                    }
                }
                OutputFormat::Alignment => {
                    // Panaroo naming: core_gene_alignment.aln
                    let path = self.output_dir.join("core_gene_alignment.aln");
                    let writer = AlignmentWriter::new();
                    writer.write_core(graph, &path)?;
                    tracing::info!("Wrote core alignment");

                    // Trim alignment with ClipKIT if requested
                    if self.trim_alignment {
                        if let Some(clipkit) = ClipKitRunner::detect() {
                            let trim_mode: TrimMode = self.trim_mode.parse().unwrap_or(TrimMode::SmartGap);
                            let trimmed_path = self.output_dir.join("core_gene_alignment.trimmed.aln");
                            match clipkit.trim(&path, &trimmed_path, trim_mode) {
                                Ok(_) => {
                                    paths.alignment = Some(trimmed_path);
                                    tracing::info!("Trimmed alignment with ClipKIT (mode: {})", trim_mode);
                                }
                                Err(e) => {
                                    tracing::warn!("ClipKIT trimming failed: {}. Using untrimmed alignment.", e);
                                }
                            }
                        } else {
                            tracing::warn!("ClipKIT not found. Install it with: pip install clipkit");
                        }
                    }

                    // Codon alignment with MACSE if requested
                    if self.codons {
                        if let Some(macse) = MacseRunner::detect() {
                            let codon_path = self.output_dir.join("core_gene_alignment.codon.aln");
                            match macse.align_codons(&path, &codon_path) {
                                Ok(_) => {
                                    paths.codon_alignment = Some(codon_path);
                                    tracing::info!("Wrote codon alignment with MACSE");
                                }
                                Err(e) => {
                                    tracing::warn!("MACSE codon alignment failed: {}. Skipping codon alignment.", e);
                                }
                            }
                        } else {
                            tracing::warn!("MACSE not found. Install it from: https://bioweb.supagro.inra.fr/macse/");
                        }
                    }

                    // BMGE filtering (if requested)
                    if self.filter_method == FilterMethod::Bmge {
                        if let Some(bmge) = BmgeRunner::detect() {
                            let bmge_path = self.output_dir.join("core_gene_alignment.BMGE.aln");
                            match bmge.filter(&path, &bmge_path, 0.2) {
                                Ok(_) => {
                                    paths.bmge_alignment = Some(bmge_path);
                                    tracing::info!("Filtered alignment with BMGE");
                                }
                                Err(e) => {
                                    tracing::warn!("BMGE filtering failed: {}", e);
                                }
                            }
                        } else {
                            tracing::warn!("BMGE not found. Install with: pip install bmge");
                        }
                    }

                    // Use trimmed path if available, otherwise original
                    if paths.alignment.is_none() {
                        paths.alignment = Some(path);
                    }
                }
                OutputFormat::Graph => {
                    // Panaroo naming: final_graph.gml
                    let path = self.output_dir.join("final_graph.gml");
                    GmlWriter::write(graph, &path)?;
                    paths.graph = Some(path);
                    tracing::info!("Wrote final GML graph");
                }
                OutputFormat::Gff => {
                    let gff_dir = self.output_dir.join("gff");
                    match write_gff_files(graph, &self.output_dir) {
                        Ok(written) => {
                            paths.gff = Some(gff_dir);
                            tracing::info!("Wrote {} GFF3 files for corrected annotations", written.len());
                        }
                        Err(e) => tracing::warn!("Failed to write GFF3 files: {}", e),
                    }
                }
                OutputFormat::Json => {
                    // Panaroo-style JSON output with gene_data.csv and pan_genome_reference.fa
                    let gene_data_path = self.output_dir.join("gene_data.csv");
                    JsonWriter::write_gene_data(graph, &graph.gene_lookup, &gene_data_path)?;
                    paths.gene_data = Some(gene_data_path);
                    tracing::info!("Wrote gene_data.csv");

                    let reference_path = self.output_dir.join("pan_genome_reference.fa");
                    JsonWriter::write_reference(graph, &reference_path)?;
                    paths.reference_fasta = Some(reference_path);
                    tracing::info!("Wrote pan_genome_reference.fa");

                    // Combined FASTA files (DNA and protein)
                    let dna_path = self.output_dir.join("combined_DNA_CDS.fasta");
                    JsonWriter::write_dna_fasta(graph, &dna_path)?;
                    paths.dna_fasta = Some(dna_path);
                    tracing::info!("Wrote combined_DNA_CDS.fasta");

                    let protein_path = self.output_dir.join("combined_protein_CDS.fasta");
                    JsonWriter::write_protein_fasta(graph, &protein_path)?;
                    paths.protein_fasta = Some(protein_path);
                    tracing::info!("Wrote combined_protein_CDS.fasta");

                    // JSON for programmatic access (different from Panaroo's JSON, but useful)
                    let json_path = self.output_dir.join("_pangenome.json");
                    JsonWriter::write_json(graph, &json_path)?;
                    paths.json = Some(json_path);
                    tracing::info!("Wrote JSON output");
                }
                OutputFormat::Parquet => {
                    #[cfg(feature = "parquet")]
                    {
                        let path = self.output_dir.join("matrix.parquet");
                        let writer = ParquetWriter::new();
                        writer.write_matrix(matrix, &path)?;
                        paths.parquet = Some(path);
                        tracing::info!("Wrote Parquet presence/absence matrix");
                    }
                    #[cfg(not(feature = "parquet"))]
                    {
                        tracing::warn!("Parquet output requires --features parquet");
                    }
                }
                OutputFormat::HtmlViz => {
                    #[cfg(feature = "viz")]
                    {
                        let path = self.output_dir.join("pangenome_viz.html");
                        let writer = HtmlVizWriter::new();
                        writer.write(graph, matrix, &path)?;
                        paths.html_viz = Some(path);
                        tracing::info!("Wrote HTML visualization");
                    }
                    #[cfg(not(feature = "viz"))]
                    {
                        tracing::warn!("HTML visualization requires --features viz");
                    }
                }
                OutputFormat::Struct => {
                    // Structural variant matrix (CSV)
                    let path = self.output_dir.join("struct_presence_absence.csv");
                    JsonWriter::write_structural_variants(graph, &path)?;
                    paths.struct_csv = Some(path);
                    tracing::info!("Wrote structural variant matrix (CSV)");
                }
                OutputFormat::SVMatrix => {
                    // Structural variant matrix (TSV)
                    let path = self.output_dir.join("struct_presence_absence.tsv");
                    let genome_names: Vec<String> = graph
                        .genomes
                        .keys()
                        .map(|g| g.as_str().to_string())
                        .collect();
                    let triplets = SVMatrixWriter::extract_triplets(graph);
                    SVMatrixWriter::new()
                        .with_genomes(genome_names)
                        .write_tsv(&triplets, &path)?;
                    paths.sv_matrix = Some(path);
                    tracing::info!("Wrote structural variant matrix (TSV)");
                }
            }
        }

        // Always write summary statistics
        let summary_path = self.output_dir.join("summary_statistics.txt");
        write_summary_stats(matrix, &summary_path, Some(graph))?;
        paths.summary_stats = Some(summary_path);
        tracing::info!("Wrote summary statistics");

        Ok(paths)
    }
}

/// Paths to generated output files.
#[derive(Debug, Clone)]
pub struct OutputPaths {
    /// Output directory
    pub output_dir: PathBuf,
    /// Gene presence/absence matrix (CSV - Roary compatible)
    pub matrix_csv: Option<PathBuf>,
    /// Gene presence/absence matrix (Rtab - Roary compatible binary)
    pub matrix_rtab: Option<PathBuf>,
    /// Gene presence/absence matrix with gene member IDs (Roary compatible)
    pub matrix_roary_csv: Option<PathBuf>,
    /// Core gene alignment
    pub alignment: Option<PathBuf>,
    /// Codon alignment (MACSE)
    pub codon_alignment: Option<PathBuf>,
    /// BMGE filtered alignment
    pub bmge_alignment: Option<PathBuf>,
    /// Final pangenome graph (GML)
    pub graph: Option<PathBuf>,
    /// GFF3 output directory (per-genome corrected annotations)
    pub gff: Option<PathBuf>,
    /// Panaroo-style reference FASTA (all genes)
    pub reference_fasta: Option<PathBuf>,
    /// Gene data CSV (links gene sequences to annotations)
    pub gene_data: Option<PathBuf>,
    /// Combined DNA CDS FASTA
    pub dna_fasta: Option<PathBuf>,
    /// Combined protein CDS FASTA
    pub protein_fasta: Option<PathBuf>,
    /// JSON output (programmatic access)
    pub json: Option<PathBuf>,
    /// JSONL streaming output
    pub jsonl: Option<PathBuf>,
    /// Structural variant matrix (CSV)
    pub struct_csv: Option<PathBuf>,
    /// Structural variant matrix (TSV)
    pub sv_matrix: Option<PathBuf>,
    /// Summary statistics file
    pub summary_stats: Option<PathBuf>,
    /// Pre-filtered graph (written before correction)
    pub pre_filt_graph: Option<PathBuf>,
    /// Parquet output files
    #[cfg(feature = "parquet")]
    pub parquet: Option<PathBuf>,
    /// HTML visualization output
    #[cfg(feature = "viz")]
    pub html_viz: Option<PathBuf>,
}
