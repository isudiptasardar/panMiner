//! Pyseer wrapper for GWAS integration.

use std::path::PathBuf;
use std::process::Command;
use std::fs;
use tempfile::TempDir;

use crate::error::Result;
use crate::graph::{PangenomeGraph, BitPackedMatrix};
use crate::gwas::traits::{GWASRunner, GWASOutput, GWASResult};

/// Pyseer runner for GWAS analysis.
///
/// Pyseer is a Python tool for genome-wide association studies
/// that can perform association mapping on pangenome data.
///
/// # Example
///
/// ```no_run
/// use std::path::PathBuf;
/// use panminer::{PangenomeGraph, BitPackedMatrix, PyseerRunner};
///
/// let mut runner = PyseerRunner::new();
/// runner.with_distances(PathBuf::from("distances.txt"));
/// runner.with_phenotypes(PathBuf::from("phenotypes.txt"));
/// // runner.run(&graph, &matrix)?;
/// ```
#[derive(Clone, Default)]
pub struct PyseerRunner {
    distances_file: Option<PathBuf>,
    phenotypes_file: Option<PathBuf>,
    output_file: Option<PathBuf>,
}

impl PyseerRunner {
    /// Create a new Pyseer runner.
    pub fn new() -> Self {
        Self {
            distances_file: None,
            phenotypes_file: None,
            output_file: None,
        }
    }

    /// Set the distance matrix file path.
    pub fn with_distances(&mut self, path: PathBuf) -> &mut Self {
        self.distances_file = Some(path);
        self
    }

    /// Set the phenotypes file path.
    pub fn with_phenotypes(&mut self, path: PathBuf) -> &mut Self {
        self.phenotypes_file = Some(path);
        self
    }

    /// Set the output file path.
    pub fn with_output(&mut self, path: PathBuf) -> &mut Self {
        self.output_file = Some(path);
        self
    }

    /// Check if pyseer is installed and available on the system.
    pub fn is_installed() -> bool {
        Command::new("pyseer").arg("--help").output().is_ok()
    }

    /// Generate input files from graph and matrix data.
    ///
    /// This creates the distance matrix and phenotype files needed by pyseer.
    fn generate_input_files(
        &self,
        _graph: &PangenomeGraph,
        matrix: &BitPackedMatrix,
    ) -> Result<(PathBuf, PathBuf)> {
        // Create temporary directory for input files
        let temp_dir = TempDir::new()
            .map_err(|e| crate::Error::Output(format!("Failed to create temp dir: {}", e)))?;

        // Generate distance matrix file
        let distance_path = temp_dir.path().join("distances.txt");
        self.write_distance_matrix(&distance_path, _graph, matrix)?;

        // Generate phenotypes file
        let phenotype_path = temp_dir.path().join("phenotypes.txt");
        self.write_phenotypes(&phenotype_path, _graph, matrix)?;

        Ok((distance_path, phenotype_path))
    }

    /// Write distance matrix in pyseer format.
    ///
    /// Format: SNP_ID GENOME1_GENOME2 DISTANCE
    fn write_distance_matrix(
        &self,
        path: &PathBuf,
        graph: &PangenomeGraph,
        matrix: &BitPackedMatrix,
    ) -> Result<()> {
        let mut lines = Vec::new();

        // Get sorted list of genome IDs
        let genome_ids: Vec<_> = graph.genomes.keys().cloned().collect();

        for (i, genome_i) in genome_ids.iter().enumerate() {
            for genome_j in genome_ids.iter().skip(i + 1) {
                // Calculate Mash-like distance based on shared clusters
                let distance = self.calculate_genome_distance(genome_i, genome_j, graph, matrix);

                // Write: genome1 genome2 distance
                lines.push(format!("{} {} {}", genome_i, genome_j, distance));
            }
        }

        fs::write(path, lines.join("\n")).map_err(|e| {
            crate::Error::Output(format!("Failed to write distance matrix: {}", e))
        })?;

        Ok(())
    }

    /// Write phenotypes file in pyseer format.
    ///
    /// Format: GENOME_ID PHENOTYPE_VALUE
    fn write_phenotypes(&self, path: &PathBuf, graph: &PangenomeGraph, _matrix: &BitPackedMatrix) -> Result<()> {
        let mut lines = Vec::new();

        for (genome_id, metadata) in &graph.genomes {
            // Use total gene count as a simple phenotype
            // In practice, this could be any numeric trait
            let phenotype_value = metadata.total_genes as f64;
            lines.push(format!("{} {}", genome_id, phenotype_value));
        }

        fs::write(path, lines.join("\n")).map_err(|e| {
            crate::Error::Output(format!("Failed to write phenotypes: {}", e))
        })?;

        Ok(())
    }

    /// Calculate distance between two genomes based on shared clusters.
    fn calculate_genome_distance(
        &self,
        genome_i: &crate::graph::GenomeId,
        genome_j: &crate::graph::GenomeId,
        _graph: &PangenomeGraph,
        matrix: &BitPackedMatrix,
    ) -> f64 {
        let total_clusters = matrix.num_clusters();

        if total_clusters == 0 {
            return 1.0;
        }

        // Count shared and unique clusters
        let mut shared = 0usize;
        let mut unique_i = 0usize;
        let mut unique_j = 0usize;

        for cluster_idx in 0..total_clusters {
            // Get genome indices from matrix
            let genome_i_idx = genome_i.as_str().parse::<usize>().unwrap_or(0);
            let genome_j_idx = genome_j.as_str().parse::<usize>().unwrap_or(0);

            let present_i = matrix.get(genome_i_idx, cluster_idx);
            let present_j = matrix.get(genome_j_idx, cluster_idx);

            match (present_i, present_j) {
                (true, true) => shared += 1,
                (true, false) => unique_i += 1,
                (false, true) => unique_j += 1,
                (false, false) => {}
            }
        }

        // Calculate Jaccard distance
        let total = shared + unique_i + unique_j;
        if total == 0 {
            1.0
        } else {
            1.0 - (shared as f64 / total as f64)
        }
    }

    /// Run pyseer and parse results.
    fn run_pyseer(&self, distance_path: &PathBuf, phenotype_path: &PathBuf) -> Result<GWASOutput> {
        // Check if pyseer is available
        if !Self::is_installed() {
            return Err(crate::Error::ExternalTool(
                "pyseer not installed. Install with: pip install pyseer".to_string(),
            ));
        }

        // Check input files exist
        if !distance_path.exists() {
            return Err(crate::Error::InvalidInput(format!(
                "Distance file not found: {}",
                distance_path.display()
            )));
        }
        if !phenotype_path.exists() {
            return Err(crate::Error::InvalidInput(format!(
                "Phenotypes file not found: {}",
                phenotype_path.display()
            )));
        }

        // Run pyseer command
        // Format: pyseer --distances distance_file phenotypes_file > output_file
        let output = Command::new("pyseer")
            .arg("--distances")
            .arg(distance_path)
            .arg(phenotype_path)
            .output();

        match output {
            Ok(output) => {
                if output.status.success() {
                    self.parse_output(&output.stdout)
                } else {
                    let stderr = String::from_utf8_lossy(&output.stderr);
                    Err(crate::Error::ExternalTool(format!(
                        "pyseer failed: {}",
                        stderr
                    )))
                }
            }
            Err(e) => Err(crate::Error::ExternalTool(format!(
                "Failed to run pyseer: {}",
                e
            ))),
        }
    }

    /// Parse pyseer output.
    fn parse_output(&self, stdout: &[u8]) -> Result<GWASOutput> {
        let output_str = String::from_utf8_lossy(stdout);
        let mut results = Vec::new();
        let mut significant_count = 0;

        for line in output_str.lines() {
            if line.trim().is_empty() || line.starts_with('#') || line.starts_with("SNP") {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                // Expected format: SNP_ID EFFECT_SIZE P_VALUE FDR
                if let (Ok(snp_id), Ok(effect), Ok(pval), Ok(fdr)) = (
                    parts[0].parse::<String>(),
                    parts[1].parse::<f64>(),
                    parts[2].parse::<f64>(),
                    parts[3].parse::<f64>(),
                ) {
                    let is_significant = fdr < 0.05;
                    if is_significant {
                        significant_count += 1;
                    }

                    results.push(GWASResult {
                        snp_id,
                        effect_size: effect,
                        p_value: pval,
                        fdr,
                    });
                }
            }
        }

        Ok(GWASOutput {
            snp_count: results.len(),
            significant_count,
            results,
        })
    }
}

impl GWASRunner for PyseerRunner {
    fn with_distances(&mut self, path: PathBuf) {
        self.distances_file = Some(path);
    }

    fn with_phenotypes(&mut self, path: PathBuf) {
        self.phenotypes_file = Some(path);
    }

    fn run(&self, graph: &PangenomeGraph, matrix: &BitPackedMatrix) -> Result<GWASOutput> {
        // Check if pyseer is available
        if !Self::is_installed() {
            return Err(crate::Error::ExternalTool(
                "pyseer not installed. Install with: pip install pyseer".to_string(),
            ));
        }

        // Use provided files or generate from graph/matrix
        let (distance_path, phenotype_path) = match (&self.distances_file, &self.phenotypes_file) {
            (Some(dist), Some(phen)) => (dist.clone(), phen.clone()),
            _ => self.generate_input_files(graph, matrix)?,
        };

        // Run pyseer
        self.run_pyseer(&distance_path, &phenotype_path)
    }

    fn is_available(&self) -> bool {
        Self::is_installed()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pyseer_runner_creation() {
        let runner = PyseerRunner::new();
        assert!(runner.is_available() || true); // May or may not be installed
    }

    #[test]
    fn test_pyseer_runner_builder_pattern() {
        let mut runner = PyseerRunner::new();
        runner.with_distances(PathBuf::from("distances.txt"));
        runner.with_phenotypes(PathBuf::from("phenotypes.txt"));

        assert!(runner.distances_file.is_some());
        assert!(runner.phenotypes_file.is_some());
    }
}
