//! MAFFT alignment runner implementation.
//!
//! Runs MAFFT subprocess for multiple sequence alignment.

use crate::clustering::alignment_traits::{AlignmentResult, AlignmentRunner, AlignmentTool};
use crate::error::{Error, Result};
use std::fmt::Write as FmtWrite;
use std::io::Write as IoWrite;
use std::process::{Command, Stdio};

/// MAFFT alignment runner.
///
/// Uses the MAFFT subprocess for multiple sequence alignment.
/// MAFFT is a fast multiple sequence alignment tool that uses FFT (Fast Fourier Transform)
/// to accelerate sequence alignment.
pub struct MafftRunner {
    /// MAFFT executable path
    executable: String,
    /// MAFFT algorithm options
    algorithm: String,
}

impl MafftRunner {
    /// Create a new MAFFT runner with default settings.
    pub fn new() -> Self {
        Self {
            executable: "mafft".to_string(),
            algorithm: "auto".to_string(), // MAFFT's auto algorithm selection
        }
    }

    /// Set the MAFFT executable path.
    pub fn with_executable(mut self, path: &str) -> Self {
        self.executable = path.to_string();
        self
    }

    /// Set the MAFFT algorithm.
    pub fn with_algorithm(mut self, algorithm: &str) -> Self {
        self.algorithm = algorithm.to_string();
        self
    }
}

impl Default for MafftRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl AlignmentRunner for MafftRunner {
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

        // Run MAFFT subprocess with FASTA input via stdin
        let mut child = Command::new(&self.executable)
            .arg("--quiet")
            .arg(self.algorithm.as_str())
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| Error::Alignment(format!("Failed to spawn MAFFT: {}", e)))?;

        // Write input to stdin
        {
            let stdin = child.stdin.as_mut().ok_or_else(|| {
                Error::Alignment("Failed to get stdin for MAFFT".to_string())
            })?;
            stdin.write_all(fasta_input.as_bytes())?;
        }

        let output = child.wait_with_output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::Alignment(format!(
                "MAFFT failed: {}",
                stderr.trim()
            )));
        }

        // Parse the aligned sequences from stdout
        let stdout = String::from_utf8_lossy(&output.stdout);
        let (aligned_sequences, alignment_length) = parse_fasta(&stdout);

        if aligned_sequences.is_empty() {
            return Err(Error::Alignment("No sequences in MAFFT output".to_string()));
        }

        Ok(AlignmentResult {
            num_sequences: aligned_sequences.len(),
            alignment_length,
            aligned_fasta: stdout.to_string(),
            tool: AlignmentTool::Mafft,
        })
    }

    fn name(&self) -> &str {
        "MAFFT Runner"
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

    // Alignment length = full concatenated sequence length (not first line length)
    let alignment_length = sequences.first().map(|(_, seq)| seq.len()).unwrap_or(0);

    (sequences, alignment_length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mafft_runner_creation() {
        let runner = MafftRunner::new();
        assert_eq!(runner.executable, "mafft");
        assert_eq!(runner.algorithm, "auto");
    }

    #[test]
    fn test_mafft_runner_custom_settings() {
        let runner = MafftRunner::new()
            .with_executable("/usr/local/bin/mafft")
            .with_algorithm("L-INS-i");
        assert_eq!(runner.executable, "/usr/local/bin/mafft");
        assert_eq!(runner.algorithm, "L-INS-i");
    }

    #[test]
    fn test_mafft_runner_default() {
        let runner = MafftRunner::default();
        assert_eq!(runner.executable, "mafft");
        assert_eq!(runner.algorithm, "auto");
    }
}
