//! Output module for PanMiner.
//!
//! Generates multiple output formats from the pangenome graph:
//! - Presence/absence matrix (CSV/TSV)
//! - Core/accessory alignments (FASTA)
//! - GML graph format
//! - JSON/JSONL
//! - Parquet (optional)
//! - Interactive HTML visualization (optional)

mod matrix;
mod alignment;
mod graph;
mod json;

pub use matrix::MatrixWriter;
pub use alignment::AlignmentWriter;
pub use graph::GmlWriter;
pub use json::JsonWriter;

use std::path::PathBuf;

use crate::config::{OutputFormat, PanminerConfig};
use crate::error::Result;
use crate::graph::PangenomeGraph;
use crate::graph::BitPackedMatrix;

/// Writes all requested output formats in parallel.
pub struct OutputWriter {
    output_dir: PathBuf,
    prefix: String,
    formats: Vec<OutputFormat>,
}

impl OutputWriter {
    /// Create a new output writer from config.
    pub fn new(config: &PanminerConfig) -> Self {
        Self {
            output_dir: config.output_dir.clone(),
            prefix: config.output_prefix.clone(),
            formats: config.outputs.iter().cloned().collect(),
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
    ) -> Result<OutputPaths> {
        self.ensure_output_dir()?;

        let mut paths = OutputPaths {
            output_dir: self.output_dir.clone(),
            matrix: None,
            alignment: None,
            graph: None,
            json: None,
        };

        // Write formats (sequential for now to collect paths)
        for format in &self.formats {
            match format {
                OutputFormat::Matrix => {
                    let path = self.output_dir.join(format!("{}_gene_presence_absence.csv", self.prefix));
                    MatrixWriter::write(matrix, &path)?;
                    paths.matrix = Some(path);
                    tracing::info!("Wrote presence/absence matrix");
                }
                OutputFormat::Alignment => {
                    let path = self.output_dir.join(format!("{}_core_alignment.fasta", self.prefix));
                    AlignmentWriter::write_core(graph, &path)?;
                    paths.alignment = Some(path);
                    tracing::info!("Wrote core alignment");
                }
                OutputFormat::Graph => {
                    let path = self.output_dir.join(format!("{}_graph.gml", self.prefix));
                    GmlWriter::write(graph, &path)?;
                    paths.graph = Some(path);
                    tracing::info!("Wrote GML graph");
                }
                OutputFormat::Json => {
                    let path = self.output_dir.join(format!("{}_pangenome.json", self.prefix));
                    JsonWriter::write(graph, matrix, &path)?;
                    paths.json = Some(path);
                    tracing::info!("Wrote JSON output");
                }
                OutputFormat::Parquet => {
                    #[cfg(feature = "parquet")]
                    {
                        tracing::info!("Parquet output: not yet implemented");
                    }
                    #[cfg(not(feature = "parquet"))]
                    {
                        tracing::warn!("Parquet output requires --features parquet");
                    }
                }
                OutputFormat::HtmlViz => {
                    #[cfg(feature = "viz")]
                    {
                        tracing::info!("HTML visualization: not yet implemented");
                    }
                    #[cfg(not(feature = "viz"))]
                    {
                        tracing::warn!("HTML visualization requires --features viz");
                    }
                }
            }
        }

        Ok(paths)
    }
}

/// Paths to generated output files.
#[derive(Debug, Clone)]
pub struct OutputPaths {
    /// Output directory
    pub output_dir: PathBuf,
    /// Presence/absence matrix
    pub matrix: Option<PathBuf>,
    /// Core alignment
    pub alignment: Option<PathBuf>,
    /// GML graph
    pub graph: Option<PathBuf>,
    /// JSON output
    pub json: Option<PathBuf>,
}
