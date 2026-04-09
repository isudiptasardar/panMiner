//! Presence/absence matrix output (CSV/TSV).

use std::path::Path;

use crate::error::Result;
use crate::graph::BitPackedMatrix;

/// Writer for presence/absence matrix output.
pub struct MatrixWriter;

impl MatrixWriter {
    /// Write the presence/absence matrix to CSV.
    pub fn write(matrix: &BitPackedMatrix, path: &Path) -> Result<()> {
        matrix.to_csv(path)
    }

    /// Write the presence/absence matrix to TSV.
    pub fn write_tsv(matrix: &BitPackedMatrix, path: &Path) -> Result<()> {
        matrix.to_tsv(path)
    }
}
