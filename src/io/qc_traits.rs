//! Quality control runner traits and definitions.
//!
//! Provides traits for running pre-processing QC tools:
//! - CheckM2: Assembly completeness and contamination scoring
//! - skani: Sparse k-mer chaining ANI for distance estimation

use crate::error::{Error, Result};
use std::path::PathBuf;

/// QC metrics for a single genome.
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct GenomeQC {
    /// Genome ID
    pub genome_id: String,
    /// Assembly completeness (0-100%)
    pub completeness: f64,
    /// Assembly contamination (0-100%)
    pub contamination: f64,
    /// Genome size in bp
    pub genome_size: u64,
    /// Number of contigs
    pub num_contigs: usize,
    /// N50 value
    pub n50: usize,
    /// Mash distance to closest reference (optional)
    pub mash_distance: Option<f64>,
    /// Whether genome passed QC
    pub passed: bool,
}

/// QC mode - how strict the QC should be.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QcMode {
    /// Strict QC - remove genomes with any issues
    Strict,
    /// Default QC - warn on issues but keep
    Default,
    /// Sensitive QC - minimal filtering
    Sensitive,
}

/// Genome distance information from ANI/distance calculations.
#[derive(Debug, Clone)]
pub struct GenomeDistance {
    /// Genome names (ordered same as matrix rows/columns)
    pub genome_names: Vec<String>,
    /// Pairwise ANI/distance matrix (0.0-1.0)
    pub distance_matrix: Vec<Vec<f64>>,
    /// Method used for distance computation
    pub method: DistanceMethod,
}

/// Method used for computing genome distances.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceMethod {
    /// skani sparse k-mer chaining (fast, robust for MAGs)
    Skani,
    /// No distance tool available
    None,
}

impl std::fmt::Display for DistanceMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DistanceMethod::Skani => write!(f, "skani"),
            DistanceMethod::None => write!(f, "none"),
        }
    }
}

impl GenomeDistance {
    /// Write the distance matrix to a CSV file.
    pub fn write_csv(&self, path: &std::path::Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        // Header
        use std::io::Write;
        write!(file, "genome")?;
        for name in &self.genome_names {
            write!(file, ",{}", name)?;
        }
        writeln!(file)?;

        // Rows
        for (i, name) in self.genome_names.iter().enumerate() {
            write!(file, "{}", name)?;
            for j in 0..self.genome_names.len() {
                write!(file, ",{:.6}", self.distance_matrix[i][j])?;
            }
            writeln!(file)?;
        }

        Ok(())
    }
}

impl QcMode {
    /// Get contamination threshold for this mode.
    pub fn contamination_threshold(&self) -> f64 {
        match self {
            QcMode::Strict => 5.0,
            QcMode::Default => 10.0,
            QcMode::Sensitive => 20.0,
        }
    }

    /// Get minimum completeness for this mode.
    pub fn min_completeness(&self) -> f64 {
        match self {
            QcMode::Strict => 90.0,
            QcMode::Default => 70.0,
            QcMode::Sensitive => 50.0,
        }
    }
}

impl Default for QcMode {
    fn default() -> Self {
        QcMode::Default
    }
}

/// Trait for QC runner implementations.
pub trait QcRunner {
    /// Run QC on the given assembly file.
    fn run_qc(&self, assembly_path: &PathBuf) -> Result<GenomeQC>;

    /// Get the name of this QC runner.
    fn name(&self) -> &str;

    /// Check if this runner is available on the system.
    fn is_available(&self) -> bool;
}

/// CheckM2-based QC runner.
///
/// CheckM2 is the successor to CheckM, providing improved completeness
/// and contamination estimation using marker sets and genome graphs.
pub struct CheckmQcRunner {
    /// Path to CheckM executable
    executable: String,
    /// Path to CheckM database
    database_path: Option<PathBuf>,
    /// QC mode for threshold selection
    mode: QcMode,
}

impl CheckmQcRunner {
    /// Create a new CheckM QC runner.
    pub fn new() -> Self {
        Self {
            executable: "checkm2".to_string(),
            database_path: None,
            mode: QcMode::Default,
        }
    }

    /// Create with custom CheckM path.
    pub fn with_path(path: impl Into<String>) -> Self {
        Self {
            executable: path.into(),
            database_path: None,
            mode: QcMode::Default,
        }
    }

    /// Create with CheckM database path.
    pub fn with_database(mut self, path: PathBuf) -> Self {
        self.database_path = Some(path);
        self
    }

    /// Set QC mode for threshold selection.
    pub fn with_mode(mut self, mode: QcMode) -> Self {
        self.mode = mode;
        self
    }

    /// Detect if CheckM2 is installed.
    pub fn detect() -> Option<Self> {
        let output = std::process::Command::new("checkm2")
            .arg("--version")
            .output();

        if output.as_ref().map(|o| o.status.success()).unwrap_or(false) {
            Some(Self::new())
        } else {
            None
        }
    }
}

impl Default for CheckmQcRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl CheckmQcRunner {
    /// Compute pairwise ANI/distance matrix for genomes.
    ///
    /// Uses skani (sparse k-mer chaining) for fast, robust distance estimation.
    /// Returns None if skani is not installed.
    pub fn compute_distance_matrix(
        &self,
        genome_paths: &[PathBuf],
    ) -> Option<GenomeDistance> {
        // Use skani (fast, robust for incomplete/MAG genomes)
        if let Some(runner) = crate::io::SkaniRunner::detect() {
            match runner.compute_ani_matrix_smart(genome_paths) {
                Ok(matrix) => {
                    let names: Vec<String> = genome_paths
                        .iter()
                        .filter_map(|p| p.file_stem().map(|s| s.to_string_lossy().to_string()))
                        .collect();
                    tracing::info!("Using skani for distance estimation");
                    return Some(GenomeDistance {
                        genome_names: names,
                        distance_matrix: matrix,
                        method: DistanceMethod::Skani,
                    });
                }
                Err(e) => {
                    tracing::warn!("skani distance computation failed: {}", e);
                }
            }
        }

        tracing::warn!("No distance estimation tool available. Install skani: conda install -c bioconda skani");
        None
    }
}

impl QcRunner for CheckmQcRunner {
    fn run_qc(&self, assembly_path: &PathBuf) -> Result<GenomeQC> {
        let temp_dir = tempfile::tempdir()?;

        // Determine command - checkm2 vs checkm
        let executable = if self.executable == "checkm2" {
            "checkm2"
        } else {
            "checkm"
        };

        let output = if executable == "checkm2" {
            // CheckM2 uses 'predict' command
            std::process::Command::new(executable)
                .arg("predict")
                .arg("-i")
                .arg(assembly_path)
                .arg("-o")
                .arg(temp_dir.path())
                .arg("-t")
                .arg("1")
                .output()
        } else {
            // CheckM (legacy) uses 'lineage' command
            std::process::Command::new(executable)
                .arg("lineage")
                .arg("-t")
                .arg("1")
                .arg(assembly_path)
                .arg(temp_dir.path())
                .output()
        }
        .map_err(|e| {
            Error::ExternalTool(format!("CheckM QC failed: {}", e))
        })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::ExternalTool(format!(
                "CheckM command failed: {}",
                stderr.trim()
            )));
        }

        // Parse CheckM output to extract QC metrics
        let stdout = String::from_utf8_lossy(&output.stdout);

        // Parse the output - format varies between checkm2 and checkm
        let (completeness, contamination, genome_size, num_contigs, n50) =
            self::parse_checkm_output(&stdout);

        Ok(GenomeQC {
            genome_id: assembly_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string(),
            completeness,
            contamination,
            genome_size,
            num_contigs,
            n50,
            mash_distance: None,
            passed: true, // Validation happens in pipeline
        })
    }

    fn name(&self) -> &str {
        "CheckM2"
    }

    fn is_available(&self) -> bool {
        // Check if checkm2 command is available
        std::process::Command::new("checkm2")
            .arg("--version")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    }
}

/// Parse CheckM output to extract QC metrics.
fn parse_checkm_output(output: &str) -> (f64, f64, u64, usize, usize) {
    let mut completeness = 100.0;
    let mut contamination = 0.0;
    let mut genome_size = 0u64;
    let mut num_contigs = 0;
    let mut n50 = 0;

    // Try to parse checkm2 output format
    for line in output.lines() {
        // CheckM2 output format: Header, then data lines with metrics
        // Looking for patterns like:
        // - Completeness: XX.XX
        // - Contamination: XX.XX
        if line.to_lowercase().contains("completeness") {
            if let Some(val) = extract_float(&line) {
                completeness = val;
            }
        }
        if line.to_lowercase().contains("contamination") {
            if let Some(val) = extract_float(&line) {
                contamination = val;
            }
        }
        if line.to_lowercase().contains("genome size") || line.to_lowercase().contains("length") {
            if let Some(val) = extract_uint(&line) {
                genome_size = val;
            }
        }
        if line.to_lowercase().contains("num contigs") || line.to_lowercase().contains("contigs") {
            if let Some(val) = extract_uint(&line) {
                num_contigs = val;
            }
        }
        if line.to_lowercase().contains("n50") {
            if let Some(val) = extract_uint(&line) {
                n50 = val;
            }
        }
    }

    // If no values found, return defaults
    (completeness, contamination, genome_size, num_contigs as usize, n50 as usize)
}

/// Extract a float value from a string.
fn extract_float(line: &str) -> Option<f64> {
    for part in line.split_whitespace() {
        if let Ok(val) = part.parse::<f64>() {
            return Some(val);
        }
    }
    None
}

/// Extract an unsigned integer from a string.
fn extract_uint(line: &str) -> Option<u64> {
    for part in line.split_whitespace() {
        if let Ok(val) = part.parse::<u64>() {
            return Some(val);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qc_mode_defaults() {
        let mode = QcMode::default();
        assert_eq!(mode, QcMode::Default);
        assert_eq!(mode.contamination_threshold(), 10.0);
    }

    #[test]
    fn test_qc_mode_strict() {
        let mode = QcMode::Strict;
        assert_eq!(mode.contamination_threshold(), 5.0);
        assert_eq!(mode.min_completeness(), 90.0);
    }

    #[test]
    fn test_qc_mode_sensitive() {
        let mode = QcMode::Sensitive;
        assert_eq!(mode.contamination_threshold(), 20.0);
        assert_eq!(mode.min_completeness(), 50.0);
    }

    #[test]
    fn test_genome_qc_default() {
        let qc = GenomeQC::default();
        assert_eq!(qc.genome_id, "");
        assert_eq!(qc.completeness, 0.0);
        assert_eq!(qc.contamination, 0.0);
        // passed defaults to false, will be set by pipeline based on thresholds
        assert!(!qc.passed);
    }
}
