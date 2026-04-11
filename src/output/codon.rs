//! Codon alignment via MACSE v2 subprocess.
//!
//! MACSE (Multiple Alignment of Coding SEquences) produces codon-aware
//! alignments that preserve reading frames. It handles frameshifts and
//! stop codons better than PAL2NAL.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::Result;

/// MACSE v2 codon alignment runner.
pub struct MacseRunner {
    /// Path to macse jar file
    jar_path: PathBuf,
    /// Path to java binary
    java_path: PathBuf,
    /// Number of threads
    threads: usize,
}

impl MacseRunner {
    /// Create a new MacseRunner with explicit paths.
    pub fn new(jar_path: PathBuf) -> Self {
        Self {
            jar_path,
            java_path: PathBuf::from("java"),
            threads: 1,
        }
    }

    /// Set the number of threads.
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Detect if MACSE is installed on the system.
    ///
    /// Looks for `macse.jar` or `macse_v2.jar` in common locations:
    /// - On PATH
    /// - In the current directory
    /// - In JAVA_HOME
    pub fn detect() -> Option<Self> {
        let java = which_java()?;

        // Try common MACSE jar locations
        let jar_candidates = [
            PathBuf::from("macse_v2.jar"),
            PathBuf::from("macse.jar"),
            PathBuf::from("/usr/share/macse/macse_v2.jar"),
            PathBuf::from("/usr/local/share/macse/macse_v2.jar"),
        ];

        for jar in &jar_candidates {
            if jar.exists() {
                return Some(Self {
                    jar_path: jar.clone(),
                    java_path: java,
                    threads: 1,
                });
            }
        }

        // Try MACSE_JAR environment variable
        if let Ok(jar_env) = std::env::var("MACSE_JAR") {
            let jar_path = PathBuf::from(jar_env);
            if jar_path.exists() {
                return Some(Self {
                    jar_path: jar_path,
                    java_path: java,
                    threads: 1,
                });
            }
        }

        None
    }

    /// Get the path to the MACSE jar file.
    pub fn path(&self) -> &Path {
        &self.jar_path
    }

    /// Align coding sequences to produce a codon-aware alignment.
    ///
    /// Uses `macse -prog alignSequences` to align nucleotide CDS.
    pub fn align_codons(&self, input_path: &Path, output_path: &Path) -> Result<PathBuf> {
        let output = Command::new(&self.java_path)
            .arg("-jar")
            .arg(&self.jar_path)
            .arg("-prog")
            .arg("alignSequences")
            .arg("-seq")
            .arg(input_path)
            .arg("-out_NT")
            .arg(output_path)
            .arg("-cpu")
            .arg(self.threads.to_string())
            .output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::Output(format!(
                "MACSE codon alignment failed: {}",
                stderr.trim()
            )));
        }

        Ok(output_path.to_path_buf())
    }

    /// Align coding sequences with a reference alignment.
    ///
    /// Uses `macse -prog enrichAlignment` to add sequences to an existing alignment.
    pub fn enrich_alignment(
        &self,
        input_alignment: &Path,
        new_sequences: &Path,
        output_path: &Path,
    ) -> Result<PathBuf> {
        let output = Command::new(&self.java_path)
            .arg("-jar")
            .arg(&self.jar_path)
            .arg("-prog")
            .arg("enrichAlignment")
            .arg("-seq")
            .arg(input_alignment)
            .arg("-seq_lr")
            .arg(new_sequences)
            .arg("-out_NT")
            .arg(output_path)
            .output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::Output(format!(
                "MACSE enrich alignment failed: {}",
                stderr.trim()
            )));
        }

        Ok(output_path.to_path_buf())
    }

    /// Get the name of this tool.
    pub fn name(&self) -> &str {
        "MACSE"
    }
}

/// Find the java binary on PATH.
fn which_java() -> Option<PathBuf> {
    which::which("java").ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_macse_runner_creation() {
        let runner = MacseRunner::new(PathBuf::from("/path/to/macse_v2.jar"));
        assert_eq!(runner.path(), Path::new("/path/to/macse_v2.jar"));
    }

    #[test]
    fn test_macse_runner_with_threads() {
        let runner = MacseRunner::new(PathBuf::from("macse.jar")).with_threads(4);
        assert_eq!(runner.threads, 4);
    }

    #[test]
    fn test_macse_name() {
        let runner = MacseRunner::new(PathBuf::from("macse.jar"));
        assert_eq!(runner.name(), "MACSE");
    }

    #[test]
    fn test_which_java() {
        // Just test that which_java doesn't panic
        let _ = which_java();
    }

    #[test]
    fn test_macse_detect() {
        // Will return None on systems without MACSE
        let result = MacseRunner::detect();
        // We don't assert it's Some because MACSE may not be installed
        // Just testing it doesn't panic
        let _ = result;
    }
}