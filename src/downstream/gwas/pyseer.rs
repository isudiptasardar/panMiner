//! Pyseer runner for downstream GWAS analysis.
//!
//! Pyseer performs gene-based GWAS with lineage structuring.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};
use crate::error::{Error, Result};

/// PyseerRunner performs gene-phenotype association testing via pyseer.
#[derive(Clone, Debug)]
pub struct PyseerRunner {
    phenotypes_file: Option<PathBuf>,
    distances_file: Option<PathBuf>,
    unitigs: bool,
    lmm: bool,
    burden: bool,
    filter_threshold: Option<u32>,
    continuous: bool,
}

impl PyseerRunner {
    /// Create a new PyseerRunner with default settings.
    pub fn new() -> Self {
        Self {
            phenotypes_file: None,
            distances_file: None,
            unitigs: false,
            lmm: false,
            burden: false,
            filter_threshold: None,
            continuous: false,
        }
    }

    /// Set the phenotypes file (TSV: genome_id<tab>phenotype).
    pub fn with_phenotypes(mut self, path: PathBuf) -> Self {
        self.phenotypes_file = Some(path);
        self
    }

    /// Enable linear mixed model (LMM) mode.
    pub fn with_lmm(mut self, enabled: bool) -> Self {
        self.lmm = enabled;
        self
    }

    /// Enable burden testing mode.
    pub fn with_burden(mut self, enabled: bool) -> Self {
        self.burden = enabled;
        self
    }

    /// Set a filter threshold for variant k-mers.
    pub fn with_filter_threshold(mut self, threshold: u32) -> Self {
        self.filter_threshold = Some(threshold);
        self
    }

    /// Enable unitigs mode.
    pub fn with_unitigs(mut self, enabled: bool) -> Self {
        self.unitigs = enabled;
        self
    }

    /// Enable continuous phenotype mode.
    pub fn with_continuous(mut self, enabled: bool) -> Self {
        self.continuous = enabled;
        self
    }

    /// Check if pyseer is installed.
    pub fn is_installed() -> bool {
        Command::new("pyseer").arg("--help").output().is_ok()
    }

    /// Compute distances from a presence/absence CSV file.
    fn compute_distances_from_pa(&self, pa_csv: &Path, output_path: &Path) -> Result<()> {
        let mut cmd = Command::new("pyseer-distance-matrices");
        cmd.arg(pa_csv);
        cmd.arg("-o").arg(output_path);

        let output = cmd.output().map_err(|e| {
            Error::ExternalTool(format!("Failed to run pyseer-distance-matrices: {}", e))
        })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::ExternalTool(format!(
                "pyseer-distance-matrices failed: {}",
                stderr
            )));
        }

        Ok(())
    }

    fn run_on_output_dir(&self, output_dir: &Path) -> Result<PyseerGWASResult> {
        let phenotypes_path = self.phenotypes_file.clone().ok_or_else(|| {
            Error::Config("Pyseer requires a phenotypes file. Use `with_phenotypes()` to set it.".to_string())
        })?;

        if !phenotypes_path.exists() {
            return Err(Error::Config(format!(
                "Phenotypes file not found: {}",
                phenotypes_path.display()
            )));
        }

        let pa_csv_path = output_dir.join("gene_presence_absence.csv");
        if !pa_csv_path.exists() {
            return Err(Error::Config(format!(
                "Gene presence/absence CSV not found at {:?}",
                pa_csv_path
            )));
        }

        let downstream_dir = output_dir.join("downstream");
        std::fs::create_dir_all(&downstream_dir)?;

        // Compute distances if not provided
        let distances_path = if let Some(ref d) = self.distances_file {
            d.clone()
        } else {
            let dist_path = downstream_dir.join("distances.gz");
            self.compute_distances_from_pa(&pa_csv_path, &dist_path)?;
            dist_path
        };

        let mut cmd = Command::new("pyseer");
        cmd.arg("--distances").arg(&distances_path);
        cmd.arg("--phenotypes").arg(&phenotypes_path);

        if self.lmm {
            cmd.arg("--lmm");
        }

        if self.burden {
            cmd.arg("--burden");
        }

        if self.continuous {
            cmd.arg("--continuous");
        }

        if let Some(threshold) = self.filter_threshold {
            cmd.arg("--filter-threshold").arg(threshold.to_string());
        }

        let output_path = downstream_dir.join("pyseer_associations.csv");
        cmd.arg("--output").arg(&output_path);

        let output = cmd.output().map_err(|e| {
            Error::ExternalTool(format!("Failed to run pyseer: {}", e))
        })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::ExternalTool(format!("pyseer failed: {}", stderr)));
        }

        let num_associations = if output_path.exists() {
            let content = std::fs::read_to_string(&output_path)?;
            content.lines().count().saturating_sub(1) // subtract header
        } else {
            0
        };

        Ok(PyseerGWASResult {
            results_path: output_path,
            num_associations,
            num_genomes: 0, // Would need to parse phenotypes to get this
        })
    }
}

impl Default for PyseerRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl DownstreamRunner for PyseerRunner {
    type Output = PyseerGWASResult;

    fn run(&self, output_dir: &Path) -> Result<Self::Output> {
        self.run_on_output_dir(output_dir)
    }

    fn name(&self) -> &str {
        "Pyseer"
    }

    fn is_available(&self) -> bool {
        Self::is_installed()
    }

    fn required_inputs(&self) -> Vec<DownstreamInput> {
        vec![
            DownstreamInput::PresenceAbsenceCsv,
            DownstreamInput::PhenotypesFile,
        ]
    }
}

/// Result of Pyseer GWAS analysis.
#[derive(Debug)]
pub struct PyseerGWASResult {
    pub results_path: PathBuf,
    pub num_associations: usize,
    pub num_genomes: usize,
}

impl DownstreamResult for PyseerGWASResult {
    fn write_to(&self, _dir: &Path) -> Result<()> {
        if !self.results_path.exists() {
            return Err(Error::InvalidInput(format!(
                "Pyseer results file not found: {}",
                self.results_path.display()
            )));
        }
        // Results are already in the downstream directory
        Ok(())
    }

    fn summary(&self) -> String {
        format!(
            "Pyseer: {} associations found, results: {}",
            self.num_associations,
            self.results_path.display()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pyseer_runner_new() {
        let runner = PyseerRunner::new();
        assert!(runner.phenotypes_file.is_none());
        assert!(!runner.lmm);
        assert!(!runner.burden);
    }

    #[test]
    fn test_pyseer_runner_builder() {
        let runner = PyseerRunner::new()
            .with_lmm(true)
            .with_burden(true)
            .with_continuous(false);
        assert!(runner.lmm);
        assert!(runner.burden);
        assert!(!runner.continuous);
    }

    #[test]
    fn test_pyseer_result_summary() {
        let result = PyseerGWASResult {
            results_path: PathBuf::from("/path/to/results.csv"),
            num_associations: 42,
            num_genomes: 100,
        };
        let summary = result.summary();
        assert!(summary.contains("42 associations"));
    }
}
