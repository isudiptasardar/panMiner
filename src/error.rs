//! Error types for PanMiner.

use std::path::PathBuf;
use thiserror::Error;

/// Main error type for PanMiner operations.
#[derive(Error, Debug)]
pub enum Error {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("GFF3 parsing error at {file}:{line}: {message}")]
    GffParse {
        file: PathBuf,
        line: usize,
        message: String,
    },

    #[error("FASTA parsing error: {0}")]
    FastaParse(String),

    #[error("Invalid sequence ID: {0}")]
    InvalidSequenceId(String),

    #[error("Cluster ID not found: {0}")]
    ClusterNotFound(String),

    #[error("Genome not found: {0}")]
    GenomeNotFound(String),

    #[error("MMseqs2 error: {0}")]
    Mmseqs(String),

    #[error("MMseqs2 not found. Install MMseqs2 for GPU-accelerated clustering.")]
    MmseqsNotFound,

    #[error("Graph construction error: {0}")]
    GraphConstruction(String),

    #[error("Output error: {0}")]
    Output(String),

    #[error("Configuration error: {0}")]
    Config(String),

    #[error("Memory error: {0}")]
    Memory(String),

    #[error("Parallel execution error: {0}")]
    Parallel(String),

    #[error("Compression error: {0}")]
    Compression(String),

    #[error("CSV error: {0}")]
    Csv(#[from] csv::Error),

    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("No genomes provided")]
    NoGenomes,

    #[error("No genes found in genome: {0}")]
    NoGenes(String),

    #[error("Feature not enabled: {0}. Compile with --features {0}")]
    FeatureNotEnabled(String),

    #[error("GPU backend error: {0}")]
    Gpu(String),

    #[error("Alignment error: {0}")]
    Alignment(String),
}

/// Result type alias for PanMiner operations.
pub type Result<T> = std::result::Result<T, Error>;

impl Error {
    /// Create a GFF parse error with context.
    pub fn gff_error(file: impl Into<PathBuf>, line: usize, message: impl Into<String>) -> Self {
        Error::GffParse {
            file: file.into(),
            line,
            message: message.into(),
        }
    }

    /// Check if this error is recoverable (can continue processing).
    pub fn is_recoverable(&self) -> bool {
        matches!(
            self,
            Error::NoGenes(_) | Error::InvalidSequenceId(_)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = Error::Mmseqs("clustering failed".to_string());
        assert!(err.to_string().contains("clustering failed"));
    }

    #[test]
    fn test_recoverable_error() {
        let err = Error::NoGenes("test.gff".to_string());
        assert!(err.is_recoverable());

        let err = Error::Io(std::io::Error::new(std::io::ErrorKind::NotFound, "file not found"));
        assert!(!err.is_recoverable());
    }
}