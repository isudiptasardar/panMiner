//! Clustal Omega alignment runner implementation.
//!
//! Runs Clustal Omega subprocess for multiple sequence alignment.

use crate::clustering::alignment_traits::{AlignmentResult, AlignmentRunner, AlignmentTool};
use crate::error::{Error, Result};
use std::fmt::Write as FmtWrite;
use std::io::Write as IoWrite;
use std::process::{Command, Stdio};

/// Clustal Omega alignment runner.
///
/// Uses the Clustal Omega subprocess for multiple sequence alignment.
/// Clustal Omega is a progressive multiple sequence alignment tool
/// that uses seeded HMM profile-profile alignments.
pub struct ClustalOmegaRunner {
    /// Clustal Omega executable path
    executable: String,
}

impl ClustalOmegaRunner {
    /// Create a new Clustal Omega runner with default settings.
    pub fn new() -> Self {
        Self {
            executable: "clustalo".to_string(),
        }
    }

    /// Set the Clustal Omega executable path.
    pub fn with_executable(mut self, path: &str) -> Self {
        self.executable = path.to_string();
        self
    }
}

impl Default for ClustalOmegaRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl AlignmentRunner for ClustalOmegaRunner {
    fn run_msa(&self, sequences: &[(String, Vec<u8>)], _tool: AlignmentTool) -> Result<AlignmentResult> {
        if sequences.is_empty() {
            return Err(Error::Alignment("No sequences provided for alignment".to_string()));
        }

        // Build FASTA input as a string
        let mut fasta_input = String::new();
        for (name, seq) in sequences {
            writeln!(&mut fasta_input, ">{}", name)
                .map_err(|e| Error::Alignment(format!("Failed to write FASTA header: {}", e)))?;
            // Write sequence in lines of 80 characters
            let seq_str = String::from_utf8_lossy(seq);
            for i in (0..seq_str.len()).step_by(80) {
                writeln!(&mut fasta_input, "{}", &seq_str[i..std::cmp::min(i + 80, seq_str.len())])
                    .map_err(|e| Error::Alignment(format!("Failed to write FASTA sequence: {}", e)))?;
            }
        }

        // Run Clustal Omega subprocess with FASTA input via stdin
        let mut child = Command::new(&self.executable)
            .arg("--quiet")
            .arg("--output-order")
            .arg("input-order")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| Error::Alignment(format!("Failed to spawn Clustal Omega: {}", e)))?;

        // Write input to stdin
        {
            let stdin = child.stdin.as_mut().ok_or_else(|| {
                Error::Alignment("Failed to get stdin for Clustal Omega".to_string())
            })?;
            stdin.write_all(fasta_input.as_bytes())?;
        }

        let output = child.wait_with_output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::Alignment(format!(
                "Clustal Omega failed: {}",
                stderr.trim()
            )));
        }

        // Parse the aligned sequences from stdout
        let stdout = String::from_utf8_lossy(&output.stdout);
        let (aligned_sequences, alignment_length) = parse_fasta(&stdout);

        if aligned_sequences.is_empty() {
            return Err(Error::Alignment("No sequences in Clustal Omega output".to_string()));
        }

        Ok(AlignmentResult {
            num_sequences: aligned_sequences.len(),
            alignment_length,
            aligned_fasta: stdout.to_string(),
            tool: AlignmentTool::ClustalOmega,
        })
    }

    fn name(&self) -> &str {
        "Clustal Omega Runner"
    }

    fn is_available(&self) -> bool {
        Command::new(&self.executable)
            .arg("--version")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    }
}

/// Parse FASTA format and return sequences with alignment length.
fn parse_fasta(fasta: &str) -> (Vec<(String, Vec<u8>)>, usize) {
    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = Vec::new();

    for line in fasta.lines() {
        if line.starts_with('>') {
            if !current_name.is_empty() {
                sequences.push((current_name, current_seq));
            }
            current_name = line[1..].to_string();
            current_seq = Vec::new();
        } else {
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                current_seq.extend_from_slice(trimmed.as_bytes());
            }
        }
    }

    if !current_name.is_empty() {
        sequences.push((current_name, current_seq));
    }

    let alignment_length = sequences.first().map(|(_, seq)| seq.len()).unwrap_or(0);

    (sequences, alignment_length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clustal_runner_creation() {
        let runner = ClustalOmegaRunner::new();
        assert_eq!(runner.executable, "clustalo");
    }

    #[test]
    fn test_clustal_runner_custom_path() {
        let runner = ClustalOmegaRunner::new().with_executable("/usr/bin/clustalo");
        assert_eq!(runner.executable, "/usr/bin/clustalo");
    }
}
