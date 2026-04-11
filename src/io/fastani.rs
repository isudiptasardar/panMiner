//! FastANI subprocess runner for pairwise ANI calculation.
//!
//! FastANI computes Average Nucleotide Identity (ANI) between genomes
//! using a fast alignment-free approach. It's used for species boundary
//! detection and genome QC.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::Result;

/// FastANI runner for computing pairwise ANI between genomes.
pub struct FastAniRunner {
    /// Path to fastani binary
    fastani_path: PathBuf,
}

impl FastAniRunner {
    /// Create a new FastAniRunner with an explicit path.
    pub fn new(fastani_path: PathBuf) -> Self {
        Self { fastani_path }
    }

    /// Detect if FastANI is installed on the system.
    ///
    /// Returns `Some(FastAniRunner)` if `fastani --version` succeeds,
    /// `None` otherwise.
    pub fn detect() -> Option<Self> {
        let path = which_fastani()?;
        Some(Self { fastani_path: path })
    }

    /// Get the path to the fastani binary.
    pub fn path(&self) -> &Path {
        &self.fastani_path
    }

    /// Compute ANI between two genomes.
    ///
    /// Returns the ANI value (0.0-1.0) or an error if computation fails.
    pub fn compute_ani(&self, query: &Path, reference: &Path) -> Result<f64> {
        let output = Command::new(&self.fastani_path)
            .arg("-q")
            .arg(query)
            .arg("-r")
            .arg(reference)
            .arg("-o")
            .arg("-") // Output to stdout
            .output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::Output(format!(
                "FastANI failed: {}",
                stderr.trim()
            )));
        }

        let stdout = String::from_utf8_lossy(&output.stdout);
        parse_fastani_ani(&stdout)
    }

    /// Compute ANI between one query genome and multiple reference genomes.
    ///
    /// Returns a list of (reference_name, ANI) pairs.
    pub fn compute_ani_one_to_many(
        &self,
        query: &Path,
        references: &[PathBuf],
    ) -> Result<Vec<(String, f64)>> {
        let ref_list = tempfile::NamedTempFile::new()?;
        let mut ref_file = std::fs::File::create(ref_list.path())?;
        for r in references {
            use std::io::Write;
            writeln!(ref_file, "{}", r.display())?;
        }
        drop(ref_file);

        let output = Command::new(&self.fastani_path)
            .arg("-q")
            .arg(query)
            .arg("--rl")
            .arg(ref_list.path())
            .arg("-o")
            .arg("-")
            .output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::Output(format!(
                "FastANI one-to-many failed: {}",
                stderr.trim()
            )));
        }

        let stdout = String::from_utf8_lossy(&output.stdout);
        parse_fastani_many(&stdout)
    }

    /// Compute an all-pairs ANI matrix.
    ///
    /// Returns a symmetric matrix where result[i][j] is the ANI between
    /// genome i and genome j.
    pub fn compute_ani_matrix(&self, genomes: &[PathBuf]) -> Result<Vec<Vec<f64>>> {
        let n = genomes.len();
        let mut matrix = vec![vec![1.0; n]; n];

        // Compute upper triangle, mirror to lower
        for i in 0..n {
            for j in (i + 1)..n {
                match self.compute_ani(&genomes[i], &genomes[j]) {
                    Ok(ani) => {
                        matrix[i][j] = ani;
                        matrix[j][i] = ani;
                    }
                    Err(e) => {
                        tracing::warn!(
                            "FastANI failed for {} vs {}: {}",
                            genomes[i].display(),
                            genomes[j].display(),
                            e
                        );
                        // Leave as 0.0 (unknown)
                        matrix[i][j] = 0.0;
                        matrix[j][i] = 0.0;
                    }
                }
            }
        }

        Ok(matrix)
    }
}

/// Find the fastani binary on PATH.
fn which_fastani() -> Option<PathBuf> {
    which::which("fastani").ok()
}

/// Parse ANI value from FastANI stdout.
///
/// FastANI output format: `query_path\tref_path\tani_value\t...\n`
fn parse_fastani_ani(output: &str) -> Result<f64> {
    for line in output.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            if let Ok(ani) = parts[2].parse::<f64>() {
                return Ok(ani / 100.0); // FastANI reports percentage, convert to 0-1
            }
        }
    }
    Err(crate::Error::Output("No ANI value found in FastANI output".to_string()))
}

/// Parse multiple ANI values from FastANI stdout.
fn parse_fastani_many(output: &str) -> Result<Vec<(String, f64)>> {
    let mut results = Vec::new();
    for line in output.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            let ref_name = PathBuf::from(parts[1])
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| parts[1].to_string());
            if let Ok(ani) = parts[2].parse::<f64>() {
                results.push((ref_name, ani / 100.0));
            }
        }
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastani_runner_creation() {
        let runner = FastAniRunner::new(PathBuf::from("/usr/bin/fastani"));
        assert_eq!(runner.path(), Path::new("/usr/bin/fastani"));
    }

    #[test]
    fn test_parse_fastani_ani() {
        let output = "/path/to/query.fna\t/path/to/ref.fna\t95.5\t1234\t1500\n";
        let ani = parse_fastani_ani(output).unwrap();
        assert!((ani - 0.955).abs() < 0.001);
    }

    #[test]
    fn test_parse_fastani_ani_no_value() {
        let output = "no\tvalid\tdata\n";
        let result = parse_fastani_ani(output);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_fastani_many() {
        let output = "/path/q1.fna\t/path/r1.fna\t97.2\t100\t200\n/path/q1.fna\t/path/r2.fna\t89.1\t50\t200\n";
        let results = parse_fastani_many(output).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].0, "r1");
        assert!((results[0].1 - 0.972).abs() < 0.001);
        assert_eq!(results[1].0, "r2");
        assert!((results[1].1 - 0.891).abs() < 0.001);
    }

    #[test]
    fn test_which_fastani() {
        // This just tests that which_fastani doesn't panic
        // It will return None on systems without fastani installed
        let _ = which_fastani();
    }
}