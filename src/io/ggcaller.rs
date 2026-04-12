//! ggCaller subprocess runner for graph-native gene calling.
//!
//! ggCaller (https://github.com/bacpop/ggCaller) is a C++ tool that performs
//! graph-based gene calling using colored de Bruijn graphs. In `--gene-finding-only`
//! mode (v1.4+), it skips Panaroo QC and outputs per-genome GFF3 annotation files
//! that can be fed directly into PanMiner's pangenome pipeline.
//!
//! This module provides a subprocess runner (`GGCallerRunner`) that detects
//! ggCaller on the system, prepares input file lists, invokes the CLI, and
//! collects the output paths.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::Result;

/// Output paths produced by a successful ggCaller run.
#[derive(Debug, Clone)]
pub struct GGCallerOutput {
    /// Directory containing per-genome GFF3 files (`<output_dir>/GFF/`).
    pub gff_dir: PathBuf,
    /// Path to the pangenome reference FASTA.
    pub gene_fasta: PathBuf,
    /// Path to the gene cluster file.
    pub cluster_file: PathBuf,
    /// Path to the gene presence/absence CSV.
    pub gene_presence_absence: PathBuf,
}

/// Subprocess runner for ggCaller gene calling.
///
/// Detects ggCaller on the system, builds input file lists, invokes the CLI
/// in `--gene-finding-only` mode, and returns structured output paths.
pub struct GGCallerRunner {
    /// Path to the ggCaller binary.
    ggcaller_path: PathBuf,
}

impl GGCallerRunner {
    /// Create a new runner with an explicit path to the ggCaller binary.
    pub fn new(ggcaller_path: PathBuf) -> Self {
        Self { ggcaller_path }
    }

    /// Detect ggCaller on the system PATH.
    ///
    /// Returns `Some(GGCallerRunner)` if the `ggcaller` binary is found,
    /// `None` otherwise.
    pub fn detect() -> Option<Self> {
        let path = which_ggcaller()?;
        Some(Self {
            ggcaller_path: path,
        })
    }

    /// Get the path to the ggCaller binary.
    pub fn path(&self) -> &Path {
        &self.ggcaller_path
    }

    /// Run ggCaller in gene-finding-only mode on a set of FASTA genomes.
    ///
    /// # Arguments
    ///
    /// * `input_fastas` — Slice of FASTA file paths to annotate.
    /// * `output_dir` — Directory where ggCaller will write its output.
    /// * `threads` — Number of threads for ggCaller (0 = auto-detect).
    ///
    /// # Returns
    ///
    /// A `GGCallerOutput` struct with paths to the GFF directory, gene FASTA,
    /// cluster file, and gene presence/absence CSV.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `input_fastas` is empty
    /// - The output directory cannot be created
    /// - The input file list cannot be written
    /// - ggCaller exits with a non-zero status
    /// - The expected GFF directory is not found after the run
    pub fn call_genes(
        &self,
        input_fastas: &[PathBuf],
        output_dir: &Path,
        threads: usize,
    ) -> Result<GGCallerOutput> {
        if input_fastas.is_empty() {
            return Err(crate::Error::InvalidInput(
                "ggCaller requires at least one input FASTA file".to_string(),
            ));
        }

        // Create output directory
        std::fs::create_dir_all(output_dir)?;

        // Write input file list (ggCaller requires a file listing FASTA paths)
        let temp_dir = tempfile::tempdir()?;
        let input_list_path = temp_dir.path().join("ggcaller_input_list.txt");
        {
            let mut list_file = std::fs::File::create(&input_list_path)?;
            for fasta_path in input_fastas {
                // Convert to absolute path for ggCaller
                let abs_path = if fasta_path.is_relative() {
                    std::env::current_dir()?.join(fasta_path)
                } else {
                    fasta_path.clone()
                };
                writeln!(list_file, "{}", abs_path.display())?;
            }
        }

        let threads_display = if threads == 0 {
            "auto".to_string()
        } else {
            threads.to_string()
        };
        tracing::info!(
            "Running ggCaller on {} genome(s) with {} thread(s)",
            input_fastas.len(),
            threads_display
        );

        // Build and run the ggCaller command
        let mut cmd = Command::new(&self.ggcaller_path);
        cmd.arg("--refs")
            .arg(&input_list_path)
            .arg("--out")
            .arg(output_dir)
            .arg("--gene-finding-only");

        if threads > 0 {
            cmd.arg("--threads").arg(threads.to_string());
        }

        let output = cmd.output().map_err(|e| {
            crate::Error::ExternalTool(format!("ggCaller failed to start: {}", e))
        })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::ExternalTool(format!(
                "ggCaller exited with error: {}",
                stderr.trim()
            )));
        }

        let gff_dir = output_dir.join("GFF");
        let gene_fasta = output_dir.join("pangenome_reference.fasta");
        let cluster_file = output_dir.join("gene_clusters.csv");
        let gene_presence_absence = output_dir.join("gene_presence_absence.csv");

        // Verify GFF directory exists
        if !gff_dir.exists() {
            return Err(crate::Error::ExternalTool(format!(
                "ggCaller output GFF directory not found: {:?}",
                gff_dir
            )));
        }

        tracing::info!("ggCaller gene calling complete. GFF output at {:?}", gff_dir);

        Ok(GGCallerOutput {
            gff_dir,
            gene_fasta,
            cluster_file,
            gene_presence_absence,
        })
    }

    /// Read the GFF directory from a ggCaller output and return sorted GFF file paths.
    ///
    /// # Arguments
    ///
    /// * `output` — A reference to a `GGCallerOutput` whose `gff_dir` will be scanned.
    ///
    /// # Returns
    ///
    /// A sorted `Vec<PathBuf>` of `.gff` files found in the GFF directory.
    pub fn parse_gff_paths(output: &GGCallerOutput) -> Result<Vec<PathBuf>> {
        if !output.gff_dir.exists() {
            return Err(crate::Error::ExternalTool(format!(
                "GFF directory does not exist: {:?}",
                output.gff_dir
            )));
        }

        let mut gff_files: Vec<PathBuf> = std::fs::read_dir(&output.gff_dir)?
            .filter_map(|entry| {
                let entry = entry.ok()?;
                let path = entry.path();
                if path.extension().and_then(|e| e.to_str()) == Some("gff") {
                    Some(path)
                } else {
                    None
                }
            })
            .collect();

        gff_files.sort();

        if gff_files.is_empty() {
            tracing::warn!("No .gff files found in {:?}", output.gff_dir);
        }

        Ok(gff_files)
    }
}

/// Find the ggCaller binary on the system PATH.
fn which_ggcaller() -> Option<PathBuf> {
    which::which("ggcaller").ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ggcaller_runner_creation() {
        let runner = GGCallerRunner::new(PathBuf::from("/usr/bin/ggcaller"));
        assert_eq!(runner.path(), Path::new("/usr/bin/ggcaller"));
    }

    #[test]
    fn test_ggcaller_detect() {
        // Just verify detect() doesn't panic — it may return None if ggcaller
        // is not installed on the test machine.
        let _ = GGCallerRunner::detect();
    }

    #[test]
    fn test_ggcaller_output_fields() {
        let output = GGCallerOutput {
            gff_dir: PathBuf::from("/tmp/ggcaller_out/GFF"),
            gene_fasta: PathBuf::from("/tmp/ggcaller_out/pangenome_reference.fasta"),
            cluster_file: PathBuf::from("/tmp/ggcaller_out/gene_clusters.csv"),
            gene_presence_absence: PathBuf::from("/tmp/ggcaller_out/gene_presence_absence.csv"),
        };

        assert!(output.gff_dir.ends_with("GFF"));
        assert!(output.gene_fasta.ends_with("pangenome_reference.fasta"));
        assert!(output.cluster_file.ends_with("gene_clusters.csv"));
        assert!(output.gene_presence_absence.ends_with("gene_presence_absence.csv"));
    }

    #[test]
    fn test_call_genes_empty_input() {
        let runner = GGCallerRunner::new(PathBuf::from("/usr/bin/ggcaller"));
        let result = runner.call_genes(&[], Path::new("/tmp/ggcaller_out"), 4);
        assert!(result.is_err());
        match result.unwrap_err() {
            crate::Error::InvalidInput(msg) => {
                assert!(msg.contains("at least one input FASTA"));
            }
            other => panic!("Expected InvalidInput, got: {:?}", other),
        }
    }

    #[test]
    fn test_which_ggcaller() {
        // Just verify which_ggcaller doesn't panic.
        let _ = which_ggcaller();
    }
}