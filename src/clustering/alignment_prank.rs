//! PRANK alignment runner implementation.
//!
//! Runs PRANK subprocess for phylogeny-aware multiple sequence alignment.

use crate::clustering::alignment_traits::{AlignmentResult, AlignmentRunner, AlignmentTool};
use crate::error::{Error, Result};
use std::fmt::Write as FmtWrite;
use std::io::Write as IoWrite;
use std::process::{Command, Stdio};

/// PRANK alignment runner.
///
/// Uses the PRANK subprocess for phylogeny-aware multiple sequence alignment.
/// PRANK is a probabilistic alignment tool that incorporates phylogenetic
/// relationships during alignment to avoid alignment artifacts.
pub struct PrankRunner {
    /// PRANK executable path
    executable: String,
    /// Whether to show phylogeny
    show_phylogeny: bool,
}

impl PrankRunner {
    /// Create a new PRANK runner with default settings.
    pub fn new() -> Self {
        Self {
            executable: "prank".to_string(),
            show_phylogeny: false,
        }
    }

    /// Set the PRANK executable path.
    pub fn with_executable(mut self, path: &str) -> Self {
        self.executable = path.to_string();
        self
    }

    /// Set whether to show phylogeny in output.
    pub fn with_show_phylogeny(mut self, show: bool) -> Self {
        self.show_phylogeny = show;
        self
    }
}

impl Default for PrankRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl AlignmentRunner for PrankRunner {
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

        // PRANK writes output to a file (not stdout). Use a temp directory
        // and -o flag to control the output path.
        let temp_dir = tempfile::tempdir()
            .map_err(|e| Error::Alignment(format!("Failed to create temp dir for PRANK: {}", e)))?;
        let output_stem = temp_dir.path().join("prank_output");

        // Build PRANK command with stdin input
        let mut cmd = Command::new(&self.executable);
        cmd.arg("-d").arg("-"); // Read from stdin
        cmd.arg("-o").arg(&output_stem); // Output file stem (PRANK adds .fas extension)
        cmd.arg("-quiet"); // Reduce output
        cmd.arg("-fasta"); // Output in FASTA format

        if self.show_phylogeny {
            cmd.arg("-showphylogeny");
        }

        // Run PRANK subprocess with FASTA input via stdin
        let mut child = cmd
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| Error::Alignment(format!("Failed to spawn PRANK: {}", e)))?;

        // Write input to stdin
        {
            let stdin = child.stdin.as_mut().ok_or_else(|| {
                Error::Alignment("Failed to get stdin for PRANK".to_string())
            })?;
            stdin.write_all(fasta_input.as_bytes())?;
        }

        let output = child.wait_with_output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::Alignment(format!(
                "PRANK failed: {}",
                stderr.trim()
            )));
        }

        // PRANK writes output to a file (default: prank_output.fas)
        let output_path = output_stem.with_extension("fas");
        let aligned_content = if output_path.exists() {
            std::fs::read_to_string(&output_path)
                .map_err(|e| Error::Alignment(format!("Failed to read PRANK output file: {}", e)))?
        } else {
            // Some PRANK versions use .fasta extension
            let alt_path = output_stem.with_extension("fasta");
            if alt_path.exists() {
                std::fs::read_to_string(&alt_path)
                    .map_err(|e| Error::Alignment(format!("Failed to read PRANK output file: {}", e)))?
            } else {
                return Err(Error::Alignment(
                    "PRANK output file not found. PRANK may have failed silently.".to_string()
                ));
            }
        };

        let (aligned_sequences, alignment_length) = parse_fasta(&aligned_content);

        if aligned_sequences.is_empty() {
            return Err(Error::Alignment("No sequences in PRANK output".to_string()));
        }

        Ok(AlignmentResult {
            num_sequences: aligned_sequences.len(),
            alignment_length,
            aligned_fasta: aligned_content,
            tool: AlignmentTool::Prank,
        })
    }

    fn name(&self) -> &str {
        "PRANK Runner"
    }

    fn is_available(&self) -> bool {
        Command::new(&self.executable)
            .arg("-version")
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
            // Push current sequence if we have one
            if !current_name.is_empty() {
                sequences.push((current_name, current_seq));
            }
            // Start new sequence
            current_name = line[1..].to_string();
            current_seq = Vec::new();
        } else {
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                current_seq.extend_from_slice(trimmed.as_bytes());
            }
        }
    }

    // Don't forget the last sequence
    if !current_name.is_empty() {
        sequences.push((current_name, current_seq));
    }

    // Alignment length = length of the first complete sequence (all lines concatenated)
    let alignment_length = sequences.first().map(|(_, seq)| seq.len()).unwrap_or(0);

    (sequences, alignment_length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prank_runner_creation() {
        let runner = PrankRunner::new();
        assert_eq!(runner.executable, "prank");
        assert!(!runner.show_phylogeny);
    }

    #[test]
    fn test_prank_runner_custom_settings() {
        let runner = PrankRunner::new()
            .with_executable("/usr/local/bin/prank")
            .with_show_phylogeny(true);
        assert_eq!(runner.executable, "/usr/local/bin/prank");
        assert!(runner.show_phylogeny);
    }
}
