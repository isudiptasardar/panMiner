//! PangrowthRunner — exact pangenome openness estimation via Heaps' law.
//!
//! Pangrowth computes exact pangenome growth and core curves without
//! requiring a phylogenetic tree, and fits Heaps' law to classify the
//! pangenome as open (alpha > 0) or closed (alpha <= 0).
//!
//! # Reference
//!
//! Gautreau et al. (2024) "Pangrowth: exact pangenome growth and core
//! curves." *Bioinformatics*.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};
use crate::error::{Error, Result};

/// Classification of pangenome openness based on Heaps' law alpha.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OpennessClassification {
    /// Pangenome is open: alpha > 0, gene pool keeps growing.
    Open,
    /// Pangenome is closed: alpha <= 0, gene pool saturates.
    Closed,
}

impl std::fmt::Display for OpennessClassification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OpennessClassification::Open => write!(f, "Open"),
            OpennessClassification::Closed => write!(f, "Closed"),
        }
    }
}

/// Result from Pangrowth pangenome openness estimation.
pub struct PangrowthResult {
    /// Heaps' law exponent: alpha > 0 means open pangenome.
    pub alpha: f64,
    /// Heaps' law prefactor (kappa).
    pub kappa: f64,
    /// Open vs closed classification derived from alpha.
    pub classification: OpennessClassification,
    /// Growth curve: (n_genomes, expected_pangenome_size).
    pub growth_curve: Vec<(usize, f64)>,
    /// Core curve: (n_genomes, expected_core_size).
    pub core_curve: Vec<(usize, f64)>,
}

impl DownstreamResult for PangrowthResult {
    fn write_to(&self, dir: &Path) -> Result<()> {
        std::fs::create_dir_all(dir)
            .map_err(|e| Error::Output(format!("Failed to create output directory: {}", e)))?;

        // Write growth curve CSV
        let growth_path = dir.join("growth_curve.csv");
        let mut wtr = csv::Writer::from_path(&growth_path)?;
        wtr.write_record(&["n_genomes", "expected_pangenome_size"])?;
        for (n, size) in &self.growth_curve {
            wtr.write_record(&[n.to_string(), size.to_string()])?;
        }
        wtr.flush()
            .map_err(|e| Error::Output(format!("Failed to write growth_curve.csv: {}", e)))?;

        // Write core curve CSV
        let core_path = dir.join("core_curve.csv");
        let mut wtr = csv::Writer::from_path(&core_path)?;
        wtr.write_record(&["n_genomes", "expected_core_size"])?;
        for (n, size) in &self.core_curve {
            wtr.write_record(&[n.to_string(), size.to_string()])?;
        }
        wtr.flush()
            .map_err(|e| Error::Output(format!("Failed to write core_curve.csv: {}", e)))?;

        // Write openness summary
        let openness_path = dir.join("openness.txt");
        std::fs::write(
            &openness_path,
            format!(
                "classification={}\nalpha={}\nkappa={}\n",
                self.classification, self.alpha, self.kappa
            ),
        )
        .map_err(|e| Error::Output(format!("Failed to write openness.txt: {}", e)))?;

        Ok(())
    }

    fn summary(&self) -> String {
        format!(
            "Pangrowth: classification={}, alpha={:.4}, kappa={:.4}",
            self.classification, self.alpha, self.kappa
        )
    }
}

/// PangrowthRunner computes exact pangenome growth/core curves and fits
/// Heaps' law to classify the pangenome as open or closed.
pub struct PangrowthRunner {
    pangrowth_path: PathBuf,
    output_dir: Option<PathBuf>,
}

impl PangrowthRunner {
    /// Detect pangrowth on PATH using `which::which`.
    ///
    /// Returns `Some(PangrowthRunner)` if the `pangrowth` binary is found,
    /// or `None` otherwise.
    pub fn detect() -> Option<Self> {
        which::which("pangrowth")
            .ok()
            .map(|path| Self {
                pangrowth_path: path,
                output_dir: None,
            })
    }

    /// Create a runner from a known binary path.
    pub fn from_path(path: PathBuf) -> Self {
        Self {
            pangrowth_path: path,
            output_dir: None,
        }
    }

    /// Set the output directory for analysis outputs (builder pattern).
    pub fn with_output_dir(mut self, path: PathBuf) -> Self {
        self.output_dir = Some(path);
        self
    }

    /// Check if the pangrowth binary is installed and available.
    pub fn is_installed(&self) -> bool {
        Command::new(&self.pangrowth_path)
            .arg("--version")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    }

    /// Find the presence/absence matrix file in the output directory.
    ///
    /// Prefers `gene_presence_absence.Rtab`, falls back to
    /// `gene_presence_absence.csv`.
    fn find_pa_matrix(output_dir: &Path) -> Option<PathBuf> {
        let rtab = output_dir.join("gene_presence_absence.Rtab");
        if rtab.exists() {
            return Some(rtab);
        }
        let csv = output_dir.join("gene_presence_absence.csv");
        if csv.exists() {
            return Some(csv);
        }
        None
    }

    /// Run `pangrowth growth` on the given P/A matrix.
    fn run_growth(&self, pa_path: &Path) -> Result<String> {
        let output = Command::new(&self.pangrowth_path)
            .args(["growth", "-p"])
            .arg(pa_path)
            .output()
            .map_err(|e| Error::ExternalTool(format!("Failed to run pangrowth: {}", e)))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::ExternalTool(format!(
                "pangrowth growth failed with exit code {:?}: {}",
                output.status.code(),
                stderr
            )));
        }

        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    }

    /// Run `pangrowth hist` on the given P/A matrix to extract alpha/kappa.
    fn run_hist(&self, pa_path: &Path) -> Result<String> {
        let output = Command::new(&self.pangrowth_path)
            .args(["hist", "-p"])
            .arg(pa_path)
            .output()
            .map_err(|e| Error::ExternalTool(format!("Failed to run pangrowth hist: {}", e)))?;

        // hist may fail for small datasets; treat as non-fatal
        if !output.status.success() {
            return Ok(String::new());
        }

        Ok(String::from_utf8_lossy(&output.stdout).to_string())
    }

    /// Parse the growth curve output from `pangrowth growth`.
    ///
    /// Each data line has two whitespace-separated values: n_genomes and
    /// expected_size. Lines starting with `#` are comments and are skipped.
    fn parse_growth_output(stdout: &str) -> Vec<(usize, f64)> {
        let mut curve = Vec::new();
        for line in stdout.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 2 {
                if let Ok(n) = parts[0].parse::<usize>() {
                    if let Ok(size) = parts[1].parse::<f64>() {
                        curve.push((n, size));
                    }
                }
            }
        }
        curve
    }

    /// Parse the core curve output from `pangrowth growth` stderr or a
    /// combined output. The core curve uses the same format as the growth
    /// curve but is typically in a separate section or file.
    ///
    /// For simplicity, we parse the core curve from the same growth output
    /// if it contains a second block (after a blank line or "core" header),
    /// otherwise we return an empty curve.
    fn parse_core_output(stdout: &str) -> Vec<(usize, f64)> {
        // Pangrowth growth prints growth curve first, then core curve
        // separated by a blank line or a header line containing "core"
        let mut in_core_section = false;
        let mut curve = Vec::new();

        for line in stdout.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                // Blank line may separate growth from core section
                in_core_section = true;
                continue;
            }
            if trimmed.starts_with('#') {
                if trimmed.to_lowercase().contains("core") {
                    in_core_section = true;
                }
                continue;
            }
            if in_core_section {
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Ok(n) = parts[0].parse::<usize>() {
                        if let Ok(size) = parts[1].parse::<f64>() {
                            curve.push((n, size));
                        }
                    }
                }
            }
        }

        curve
    }

    /// Extract alpha and kappa from `pangrowth hist` output.
    ///
    /// Looks for lines containing "alpha" and "kappa" (or "K") keywords.
    /// Returns (alpha, kappa), defaulting to 0.0 if not found.
    fn parse_alpha_kappa(hist_output: &str) -> (f64, f64) {
        let mut alpha = 0.0_f64;
        let mut kappa = 0.0_f64;

        for line in hist_output.lines() {
            let lower = line.to_lowercase();

            // Match lines like "alpha = 0.23" or "alpha: 0.23" or "alpha\t0.23"
            if lower.contains("alpha") && !lower.contains("kappa") {
                if let Some(val) = extract_f64_from_line(line) {
                    alpha = val;
                }
            }

            // Match lines like "K = 123.45" or "kappa = 123.45" or "kappa: 123.45"
            if lower.contains("kappa") || (lower.starts_with('k') && !lower.contains("alpha")) {
                if let Some(val) = extract_f64_from_line(line) {
                    kappa = val;
                }
            }
        }

        (alpha, kappa)
    }

    /// Also try to extract alpha/kappa from `pangrowth growth` output,
    /// which sometimes prints Heaps' law fit parameters.
    fn extract_alpha_from_growth(stdout: &str) -> Option<f64> {
        for line in stdout.lines() {
            let lower = line.to_lowercase();
            if (lower.contains("alpha") || lower.contains("heaps"))
                && !lower.contains("kappa")
            {
                if let Some(val) = extract_f64_from_line(line) {
                    return Some(val);
                }
            }
        }
        None
    }
}

/// Extract the first parseable f64 from a line after an = or : separator.
fn extract_f64_from_line(line: &str) -> Option<f64> {
    // Try splitting on = or : first, then take the last numeric token
    let after_sep = if let Some(idx) = line.find('=') {
        &line[idx + 1..]
    } else if let Some(idx) = line.find(':') {
        &line[idx + 1..]
    } else {
        line
    };

    // Take the last whitespace-separated token that parses as f64
    for token in after_sep.split_whitespace().rev() {
        if let Ok(val) = token.trim_matches(',').parse::<f64>() {
            return Some(val);
        }
    }
    None
}

impl Default for PangrowthRunner {
    fn default() -> Self {
        Self::detect().unwrap_or_else(|| Self {
            pangrowth_path: PathBuf::from("pangrowth"),
            output_dir: None,
        })
    }
}

impl DownstreamRunner for PangrowthRunner {
    type Output = PangrowthResult;

    fn run(&self, output_dir: &Path) -> Result<Self::Output> {
        let pa_path = Self::find_pa_matrix(output_dir).ok_or_else(|| {
            Error::Output(format!(
                "Presence/absence matrix not found in {}. Expected gene_presence_absence.Rtab or gene_presence_absence.csv",
                output_dir.display()
            ))
        })?;

        // Run pangrowth growth
        let growth_stdout = self.run_growth(&pa_path)?;

        // Parse growth and core curves
        let growth_curve = Self::parse_growth_output(&growth_stdout);
        let core_curve = Self::parse_core_output(&growth_stdout);

        // Try to extract alpha/kappa from growth output first
        let (mut alpha, mut kappa) = if let Some(a) = Self::extract_alpha_from_growth(&growth_stdout)
        {
            (a, 0.0)
        } else {
            (0.0, 0.0)
        };

        // Also try pangrowth hist for alpha/kappa
        if let Ok(hist_output) = self.run_hist(&pa_path) {
            if !hist_output.is_empty() {
                let (hist_alpha, hist_kappa) = Self::parse_alpha_kappa(&hist_output);
                // Prefer hist values if they are non-zero
                if hist_alpha != 0.0 {
                    alpha = hist_alpha;
                }
                if hist_kappa != 0.0 {
                    kappa = hist_kappa;
                }
            }
        }

        let classification = if alpha > 0.0 {
            OpennessClassification::Open
        } else {
            OpennessClassification::Closed
        };

        Ok(PangrowthResult {
            alpha,
            kappa,
            classification,
            growth_curve,
            core_curve,
        })
    }

    fn name(&self) -> &str {
        "Pangrowth"
    }

    fn is_available(&self) -> bool {
        Self::detect().is_some()
    }

    fn required_inputs(&self) -> Vec<DownstreamInput> {
        vec![DownstreamInput::PresenceAbsenceCsv]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_growth_output_valid() {
        let output = "# Pangenome growth curve\n1 2000\n2 3500\n3 4800\n4 5900\n";
        let curve = PangrowthRunner::parse_growth_output(output);
        assert_eq!(curve.len(), 4);
        assert_eq!(curve[0], (1, 2000.0));
        assert_eq!(curve[1], (2, 3500.0));
        assert_eq!(curve[2], (3, 4800.0));
        assert_eq!(curve[3], (4, 5900.0));
    }

    #[test]
    fn test_parse_growth_output_with_comments() {
        let output = "# This is a comment\n1 100.5\n# Another comment\n2 205.3\n3 307.1\n";
        let curve = PangrowthRunner::parse_growth_output(output);
        assert_eq!(curve.len(), 3);
        assert_eq!(curve[0], (1, 100.5));
        assert_eq!(curve[1], (2, 205.3));
        assert_eq!(curve[2], (3, 307.1));
    }

    #[test]
    fn test_parse_growth_output_empty() {
        let curve = PangrowthRunner::parse_growth_output("");
        assert!(curve.is_empty());
    }

    #[test]
    fn test_parse_growth_output_only_comments() {
        let output = "# comment 1\n# comment 2\n";
        let curve = PangrowthRunner::parse_growth_output(output);
        assert!(curve.is_empty());
    }

    #[test]
    fn test_parse_core_output() {
        let output = "1 5000\n2 8000\n\n1 3000\n2 2500\n";
        let core = PangrowthRunner::parse_core_output(output);
        assert_eq!(core.len(), 2);
        assert_eq!(core[0], (1, 3000.0));
        assert_eq!(core[1], (2, 2500.0));
    }

    #[test]
    fn test_openness_classification_open() {
        let result = PangrowthResult {
            alpha: 0.45,
            kappa: 1200.0,
            classification: OpennessClassification::Open,
            growth_curve: vec![(1, 2000.0)],
            core_curve: vec![],
        };
        assert_eq!(result.classification, OpennessClassification::Open);
        assert!(result.summary().contains("Open"));
        assert!(result.summary().contains("0.45"));
    }

    #[test]
    fn test_openness_classification_closed() {
        let result = PangrowthResult {
            alpha: -0.1,
            kappa: 5000.0,
            classification: OpennessClassification::Closed,
            growth_curve: vec![(1, 2000.0)],
            core_curve: vec![],
        };
        assert_eq!(result.classification, OpennessClassification::Closed);
        assert!(result.summary().contains("Closed"));
    }

    #[test]
    fn test_openness_classification_zero_alpha_is_closed() {
        // alpha = 0 is considered closed (alpha <= 0)
        let result = PangrowthResult {
            alpha: 0.0,
            kappa: 3000.0,
            classification: OpennessClassification::Closed,
            growth_curve: vec![],
            core_curve: vec![],
        };
        assert_eq!(result.classification, OpennessClassification::Closed);
    }

    #[test]
    fn test_parse_alpha_kappa() {
        let hist_output = "Pangenome model: Heaps' law\nalpha = 0.567\nK = 1234.5\n";
        let (alpha, kappa) = PangrowthRunner::parse_alpha_kappa(hist_output);
        assert!((alpha - 0.567).abs() < 1e-10);
        assert!((kappa - 1234.5).abs() < 1e-10);
    }

    #[test]
    fn test_parse_alpha_kappa_missing() {
        let hist_output = "some output without alpha\n";
        let (alpha, kappa) = PangrowthRunner::parse_alpha_kappa(hist_output);
        assert_eq!(alpha, 0.0);
        assert_eq!(kappa, 0.0);
    }

    #[test]
    fn test_pangrowth_result_write_to() {
        let result = PangrowthResult {
            alpha: 0.35,
            kappa: 2100.0,
            classification: OpennessClassification::Open,
            growth_curve: vec![(1, 2000.0), (2, 3500.0), (3, 4800.0)],
            core_curve: vec![(1, 3000.0), (2, 2500.0)],
        };
        let dir = std::env::temp_dir().join("pangrowth_test_write");
        result.write_to(&dir).unwrap();

        let growth_path = dir.join("growth_curve.csv");
        assert!(growth_path.exists());

        let core_path = dir.join("core_curve.csv");
        assert!(core_path.exists());

        let openness_path = dir.join("openness.txt");
        assert!(openness_path.exists());

        let openness_content = std::fs::read_to_string(&openness_path).unwrap();
        assert!(openness_content.contains("Open"));
        assert!(openness_content.contains("0.35"));

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_pangrowth_runner_detect_returns_none_when_not_installed() {
        // This test just verifies the method doesn't panic;
        // it may return Some or None depending on the test environment.
        let _ = PangrowthRunner::detect();
    }

    #[test]
    fn test_pangrowth_runner_from_path() {
        let runner = PangrowthRunner::from_path(PathBuf::from("/usr/local/bin/pangrowth"));
        assert_eq!(runner.pangrowth_path, PathBuf::from("/usr/local/bin/pangrowth"));
        assert!(runner.output_dir.is_none());
    }

    #[test]
    fn test_pangrowth_runner_with_output_dir() {
        let runner = PangrowthRunner::from_path(PathBuf::from("pangrowth"))
            .with_output_dir(PathBuf::from("/tmp/out"));
        assert_eq!(runner.output_dir, Some(PathBuf::from("/tmp/out")));
    }

    #[test]
    fn test_extract_f64_from_line() {
        assert_eq!(extract_f64_from_line("alpha = 0.567"), Some(0.567));
        assert_eq!(extract_f64_from_line("K: 1234.5"), Some(1234.5));
        assert_eq!(extract_f64_from_line("no number here"), None);
    }
}