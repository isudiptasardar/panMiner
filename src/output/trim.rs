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

/// BMGE (Block Mapping and Gathering with Entropy) alignment filter runner.
///
/// BMGE filters poorly aligned columns from MSAs using entropy-based scoring.
/// It runs via Python/Biopython.
///
/// Reference: Criscuolo & Gribaldo, BMC Evolutionary Biology 10, 210 (2010).
pub struct BmgeRunner {
    python_path: PathBuf,
}

impl BmgeRunner {
    /// Create a new BmgeRunner with an explicit Python path.
    pub fn new(python_path: PathBuf) -> Self {
        Self { python_path }
    }

    /// Detect if BMGE is available via Python/Biopython.
    pub fn detect() -> Option<Self> {
        let python = if which::which("python3").is_ok() {
            PathBuf::from("python3")
        } else if which::which("python").is_ok() {
            PathBuf::from("python")
        } else {
            return None;
        };

        let output = std::process::Command::new(&python)
            .arg("-c")
            .arg("import bmge")
            .output()
            .ok()?;

        if output.status.success() {
            Some(Self { python_path: python })
        } else {
            None
        }
    }

    /// Get the runner name.
    pub fn name(&self) -> &str {
        "BMGE"
    }

    /// Filter an alignment using BMGE.
    pub fn filter(
        &self,
        input_path: &Path,
        output_path: &Path,
        gap_threshold: f64,
    ) -> crate::error::Result<PathBuf> {
        let script = format!(
            "import sys\nfrom Bio import AlignIO\nfrom bmge import bmge as bmge_filter\nalignment = AlignIO.read(sys.argv[1], 'fasta')\nfiltered = bmge_filter(alignment, gap_threshold={})\nAlignIO.write(filtered, sys.argv[2], 'fasta')\n",
            gap_threshold
        );

        let output = std::process::Command::new(&self.python_path)
            .arg("-c")
            .arg(&script)
            .arg(input_path)
            .arg(output_path)
            .output()
            .map_err(|e| crate::Error::Output(format!("BMGE filter failed: {}", e)))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::Output(format!(
                "BMGE filtering failed: {}. Install with: pip install bmge",
                stderr.trim()
            )));
        }

        Ok(output_path.to_path_buf())
    }
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

    #[test]
    fn test_bmge_runner_creation() {
        let runner = BmgeRunner::new(PathBuf::from("/usr/bin/python3"));
        assert_eq!(runner.name(), "BMGE");
    }

    #[test]
    fn test_bmge_detect() {
        let _ = BmgeRunner::detect();
    }
}