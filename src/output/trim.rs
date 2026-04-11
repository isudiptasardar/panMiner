//! Alignment trimming via ClipKIT subprocess.
//!
//! ClipKIT is a multiple sequence alignment trimming tool that removes
//! poorly aligned positions. It's superior to BMGE per the 138K alignment
//! benchmark (https://gitlab.com/LBC-SciBio/clipkit).

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::Result;

/// Trimming mode for ClipKIT.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrimMode {
    /// Smart gap trimming (default, recommended)
    SmartGap,
    /// Gappyout trimming (aggressive)
    Gappyout,
    /// Strict trimming (conservative)
    Strict,
}

impl std::fmt::Display for TrimMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TrimMode::SmartGap => write!(f, "smart-gap"),
            TrimMode::Gappyout => write!(f, "gappyout"),
            TrimMode::Strict => write!(f, "strict"),
        }
    }
}

impl std::str::FromStr for TrimMode {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "smart-gap" | "smart_gap" | "smart" => Ok(TrimMode::SmartGap),
            "gappyout" | "gappy" => Ok(TrimMode::Gappyout),
            "strict" => Ok(TrimMode::Strict),
            _ => Err(format!("Unknown trim mode: {}. Use smart-gap, gappyout, or strict.", s)),
        }
    }
}

/// ClipKIT alignment trimming runner.
pub struct ClipKitRunner {
    /// Path to clipkit binary
    clipkit_path: PathBuf,
}

impl ClipKitRunner {
    /// Create a new ClipKitRunner with an explicit path.
    pub fn new(clipkit_path: PathBuf) -> Self {
        Self { clipkit_path }
    }

    /// Detect if ClipKIT is installed on the system.
    pub fn detect() -> Option<Self> {
        let path = which_clipkit()?;
        Some(Self { clipkit_path: path })
    }

    /// Get the path to the clipkit binary.
    pub fn path(&self) -> &Path {
        &self.clipkit_path
    }

    /// Trim an alignment file using ClipKIT.
    ///
    /// Returns the path to the trimmed output file.
    pub fn trim(
        &self,
        input_path: &Path,
        output_path: &Path,
        mode: TrimMode,
    ) -> Result<PathBuf> {
        let mode_flag = match mode {
            TrimMode::SmartGap => "smart-gap",
            TrimMode::Gappyout => "gappyout",
            TrimMode::Strict => "strict",
        };

        let output = Command::new(&self.clipkit_path)
            .arg(input_path)
            .arg("-m")
            .arg(mode_flag)
            .arg("-o")
            .arg(output_path)
            .output()?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::Output(format!(
                "ClipKIT trimming failed: {}",
                stderr.trim()
            )));
        }

        Ok(output_path.to_path_buf())
    }

    /// Get the name of this tool.
    pub fn name(&self) -> &str {
        "ClipKIT"
    }
}

/// Find the clipkit binary on PATH.
fn which_clipkit() -> Option<PathBuf> {
    which::which("clipkit").ok().or_else(|| which::which("clipkit.py").ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clipkit_runner_creation() {
        let runner = ClipKitRunner::new(PathBuf::from("/usr/bin/clipkit"));
        assert_eq!(runner.path(), Path::new("/usr/bin/clipkit"));
    }

    #[test]
    fn test_trim_mode_display() {
        assert_eq!(TrimMode::SmartGap.to_string(), "smart-gap");
        assert_eq!(TrimMode::Gappyout.to_string(), "gappyout");
        assert_eq!(TrimMode::Strict.to_string(), "strict");
    }

    #[test]
    fn test_trim_mode_from_str() {
        assert_eq!("smart-gap".parse::<TrimMode>().unwrap(), TrimMode::SmartGap);
        assert_eq!("gappyout".parse::<TrimMode>().unwrap(), TrimMode::Gappyout);
        assert_eq!("strict".parse::<TrimMode>().unwrap(), TrimMode::Strict);
        assert!("invalid".parse::<TrimMode>().is_err());
    }

    #[test]
    fn test_which_clipkit() {
        // Just test that which_clipkit doesn't panic
        let _ = which_clipkit();
    }

    #[test]
    fn test_clipkit_name() {
        let runner = ClipKitRunner::new(PathBuf::from("clipkit"));
        assert_eq!(runner.name(), "ClipKIT");
    }
}