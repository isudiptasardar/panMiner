//! PanstripeRunner — phylogenetically-informed gene gain/loss rate estimation.
//!
//! Panstripe uses generalized linear models to estimate gene gain and loss rates
//! while controlling for population structure and sampling bias.
//!
//! # Reference
//!
//! Tonkin-Hill et al. (2023) "Panstripe: phylogenetically-informed gene
//! gain and loss rates with lineage-specific detection." *Genome Research*.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};
use crate::error::{Error, Result};

/// PanstripeRunner estimates gene gain and loss rates using phylogenetic GLMs.
pub struct PanstripeRunner {
    tree_file: Option<PathBuf>,
    output_dir: Option<PathBuf>,
}

impl PanstripeRunner {
    /// Create a new PanstripeRunner with default settings.
    pub fn new() -> Self {
        Self {
            tree_file: None,
            output_dir: None,
        }
    }

    /// Set the phylogenetic tree file (Newick format).
    pub fn with_tree(mut self, path: PathBuf) -> Self {
        self.tree_file = Some(path);
        self
    }

    /// Set the output directory for analysis outputs.
    #[allow(dead_code)]
    pub fn with_output_dir(mut self, path: PathBuf) -> Self {
        self.output_dir = Some(path);
        self
    }

    /// Detect if Rscript + panstripe R package are available.
    pub fn detect() -> Option<Self> {
        if which::which("Rscript").is_err() {
            return None;
        }
        if !Self::check_panstripe_r_package() {
            return None;
        }
        Some(Self::new())
    }

    /// Check if the panstripe R package is installed.
    fn check_panstripe_r_package() -> bool {
        let output = Command::new("Rscript")
            .args(["-e", "library(panstripe)"])
            .output();
        match output {
            Ok(out) => out.status.success(),
            Err(_) => false,
        }
    }

    /// Build the R script for running panstripe analysis.
    fn build_r_script(pa_path: &Path, tree_path: &Path) -> String {
        format!(
            r#"
library(panstripe)
pa <- read.delim("{}", row.names=1)
tree <- read.tree("{}")
result <- panstripe(pa, tree)
write.csv(coef(result), "panstripe_rates.csv")
sink("panstripe_summary.txt")
print(summary(result))
sink()
writeLines(as.character(result@convergence), "panstripe_convergence.txt")
cat("Panstripe analysis complete.\n")
"#,
            pa_path.display(),
            tree_path.display()
        )
    }

    fn execute_r_script(&self, work_dir: &Path, pa_path: &Path, tree_path: &Path) -> Result<()> {
        let r_script = Self::build_r_script(pa_path, tree_path);
        let script_path = work_dir.join("panstripe_run.R");
        std::fs::write(&script_path, &r_script)
            .map_err(|e| Error::Output(format!("Failed to write R script: {}", e)))?;
        let output = Command::new("Rscript")
            .arg(&script_path)
            .current_dir(work_dir)
            .output()
            .map_err(|e| Error::ExternalTool(format!("Failed to run Rscript: {}", e)))?;
        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::ExternalTool(format!(
                "Rscript failed with exit code {:?}: {}",
                output.status.code(),
                stderr
            )));
        }
        Ok(())
    }

    fn parse_rates_csv(path: &Path) -> Result<PanstripeResult> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| Error::Output(format!("Failed to read panstripe_rates.csv: {}", e)))?;
        let mut gain_rate = 0.0;
        let mut loss_rate = 0.0;
        let mut alpha = 0.0;
        let mut convergence = 0.0;
        let mut r_squared = 0.0;
        for line in content.lines() {
            if line.starts_with(',') || line.starts_with('"') {
                continue;
            }
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() < 2 {
                continue;
            }
            let row_name = parts[0].trim_matches('"').trim();
            if let Ok(value) = parts[1].trim().parse::<f64>() {
                match row_name {
                    "gain" | "gain_rate" => gain_rate = value,
                    "loss" | "loss_rate" => loss_rate = value,
                    "alpha" => alpha = value,
                    "convergence" | "converged" => convergence = value,
                    "r_squared" | "R2" | "adj.r.squared" => r_squared = value,
                    _ => {}
                }
            }
        }
        let conv_path = path.with_file_name("panstripe_convergence.txt");
        if conv_path.exists() {
            if let Ok(conv_str) = std::fs::read_to_string(&conv_path) {
                if let Ok(conv_val) = conv_str.trim().parse::<f64>() {
                    convergence = conv_val;
                }
            }
        }
        Ok(PanstripeResult { gain_rate, loss_rate, alpha, convergence, r_squared })
    }
}

impl Default for PanstripeRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl DownstreamRunner for PanstripeRunner {
    type Output = PanstripeResult;

    fn run(&self, output_dir: &Path) -> Result<Self::Output> {
        let tree_path = self.tree_file.clone().unwrap_or_else(|| output_dir.join("tree.nwk"));
        if !tree_path.exists() {
            return Err(Error::Output(format!("Phylogenetic tree file not found: {}", tree_path.display())));
        }
        let pa_path = output_dir.join("gene_presence_absence.Rtab");
        let pa_csv_path = output_dir.join("gene_presence_absence.csv");
        let pa_file = if pa_path.exists() {
            pa_path.as_path()
        } else if pa_csv_path.exists() {
            pa_csv_path.as_path()
        } else {
            return Err(Error::Output(format!(
                "Presence/absence matrix not found in {}. Expected gene_presence_absence.Rtab or gene_presence_absence.csv",
                output_dir.display()
            )));
        };
        let temp_dir = tempfile::tempdir()
            .map_err(|e| Error::Output(format!("Failed to create temp directory: {}", e)))?;
        let work_dir = temp_dir.path();
        let temp_pa = work_dir.join("gene_presence_absence.Rtab");
        let temp_tree = work_dir.join("tree.nwk");
        std::fs::copy(pa_file, &temp_pa)
            .map_err(|e| Error::Output(format!("Failed to copy P/A matrix: {}", e)))?;
        std::fs::copy(&tree_path, &temp_tree)
            .map_err(|e| Error::Output(format!("Failed to copy tree file: {}", e)))?;
        self.execute_r_script(work_dir, &temp_pa, &temp_tree)?;
        let rates_path = work_dir.join("panstripe_rates.csv");
        if !rates_path.exists() {
            return Err(Error::Output("panstripe_rates.csv not found after R script execution".to_string()));
        }
        let result = Self::parse_rates_csv(&rates_path)?;
        let downstream_dir = output_dir.join("downstream");
        std::fs::create_dir_all(&downstream_dir)
            .map_err(|e| Error::Output(format!("Failed to create downstream directory: {}", e)))?;
        let dest_rates = downstream_dir.join("panstripe_rates.csv");
        std::fs::copy(&rates_path, &dest_rates)
            .map_err(|e| Error::Output(format!("Failed to copy rates CSV: {}", e)))?;
        let summary_path = work_dir.join("panstripe_summary.txt");
        if summary_path.exists() {
            let _ = std::fs::copy(&summary_path, downstream_dir.join("panstripe_summary.txt"));
        }
        let conv_path = work_dir.join("panstripe_convergence.txt");
        if conv_path.exists() {
            let _ = std::fs::copy(&conv_path, downstream_dir.join("panstripe_convergence.txt"));
        }
        Ok(result)
    }

    fn name(&self) -> &str { "Panstripe" }

    fn is_available(&self) -> bool { Self::detect().is_some() }

    fn required_inputs(&self) -> Vec<DownstreamInput> {
        vec![DownstreamInput::PresenceAbsenceCsv, DownstreamInput::PhylogeneticTree]
    }
}

pub struct PanstripeResult {
    pub gain_rate: f64,
    pub loss_rate: f64,
    pub alpha: f64,
    pub convergence: f64,
    pub r_squared: f64,
}

impl DownstreamResult for PanstripeResult {
    fn write_to(&self, dir: &Path) -> Result<()> {
        std::fs::create_dir_all(dir)
            .map_err(|e| Error::Output(format!("Failed to create output directory: {}", e)))?;
        let csv_path = dir.join("panstripe_rates.csv");
        let mut wtr = csv::Writer::from_path(&csv_path)?;
        wtr.write_record(&["parameter", "value"])?;
        wtr.write_record(&["gain_rate", &self.gain_rate.to_string()])?;
        wtr.write_record(&["loss_rate", &self.loss_rate.to_string()])?;
        wtr.write_record(&["alpha", &self.alpha.to_string()])?;
        wtr.write_record(&["convergence", &self.convergence.to_string()])?;
        wtr.write_record(&["r_squared", &self.r_squared.to_string()])?;
        wtr.flush()
            .map_err(|e| Error::Output(format!("Failed to write panstripe_rates.csv: {}", e)))?;
        Ok(())
    }

    fn summary(&self) -> String {
        format!(
            "Panstripe: gain_rate={}, loss_rate={}, alpha={}, convergence={}, r_squared={}",
            self.gain_rate, self.loss_rate, self.alpha, self.convergence, self.r_squared
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_panstripe_runner_new() {
        let runner = PanstripeRunner::new();
        assert!(runner.tree_file.is_none());
        assert!(runner.output_dir.is_none());
    }

    #[test]
    fn test_panstripe_runner_with_tree() {
        let runner = PanstripeRunner::new().with_tree(PathBuf::from("/path/to/tree.nwk"));
        assert!(runner.tree_file.is_some());
        assert_eq!(runner.tree_file.unwrap(), PathBuf::from("/path/to/tree.nwk"));
    }

    #[test]
    fn test_panstripe_result_summary() {
        let result = PanstripeResult { gain_rate: 0.05, loss_rate: 0.03, alpha: 1.67, convergence: 1.0, r_squared: 0.85 };
        let summary = result.summary();
        assert!(summary.contains("gain_rate=0.05"));
        assert!(summary.contains("loss_rate=0.03"));
        assert!(summary.contains("alpha=1.67"));
    }

    #[test]
    fn test_build_r_script() {
        let pa_path = Path::new("/data/pa.Rtab");
        let tree_path = Path::new("/data/tree.nwk");
        let script = PanstripeRunner::build_r_script(pa_path, tree_path);
        assert!(script.contains("library(panstripe)"));
        assert!(script.contains(r#"read.delim("/data/pa.Rtab""#));
        assert!(script.contains(r#"read.tree("/data/tree.nwk")"#));
        assert!(script.contains("panstripe(pa, tree)"));
    }

    #[test]
    fn test_downstream_result_write() {
        let result = PanstripeResult { gain_rate: 0.1, loss_rate: 0.05, alpha: 2.0, convergence: 1.0, r_squared: 0.9 };
        let dir = std::env::temp_dir().join("panstripe_test");
        result.write_to(&dir).unwrap();
        let csv_path = dir.join("panstripe_rates.csv");
        assert!(csv_path.exists());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
