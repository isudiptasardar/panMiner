//! GWAS runner trait definition.

use std::path::PathBuf;
use crate::error::Result;
use crate::graph::PangenomeGraph;

/// Trait for GWAS (Genome-Wide Association Studies) runners.
///
/// Implementations can use different tools (Pyseer, etc.)
/// to perform association studies on pangenome data.
pub trait GWASRunner {
    /// Set the distance matrix file path.
    fn with_distances(&mut self, path: PathBuf);

    /// Set the phenotypes file path.
    fn with_phenotypes(&mut self, path: PathBuf);

    /// Run GWAS analysis on the given graph and matrix.
    ///
    /// # Arguments
    ///
    /// * `graph` - The pangenome graph
    /// * `matrix` - The bit-packed presence/absence matrix
    ///
    /// # Returns
    ///
    /// The GWAS results including SNP IDs, effect sizes, p-values, and FDR.
    fn run(&self, graph: &PangenomeGraph, matrix: &crate::graph::BitPackedMatrix) -> Result<GWASOutput>;

    /// Check if this GWAS runner is available on the system.
    fn is_available(&self) -> bool;
}

/// Output from a GWAS analysis.
#[derive(Debug, Clone, Default)]
pub struct GWASOutput {
    /// Number of SNPs tested
    pub snp_count: usize,
    /// Number of significant SNPs (FDR < 0.05)
    pub significant_count: usize,
    /// List of GWAS results
    pub results: Vec<GWASResult>,
}

/// Single GWAS result for a SNP.
#[derive(Debug, Clone, Default)]
pub struct GWASResult {
    /// SNP/cluster identifier
    pub snp_id: String,
    /// Effect size (beta coefficient)
    pub effect_size: f64,
    /// P-value
    pub p_value: f64,
    /// False Discovery Rate (FDR)
    pub fdr: f64,
}
