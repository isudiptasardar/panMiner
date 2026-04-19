//! Clustering trait definition.

use crate::error::Result;
use crate::graph::{Gene, GeneCluster};

/// Trait for gene clustering implementations.
///
/// Implementations can use different backends (MMseqs2-GPU, CPU, etc.)
/// but all produce the same output format.
pub trait Clusterer {
    /// Cluster genes by sequence similarity.
    ///
    /// # Arguments
    ///
    /// * `genes` - Genes to cluster
    /// * `identity_threshold` - Minimum identity for clustering (0.0-1.0)
    /// * `len_dif_percent` - Minimum length coverage for clustering (0.0-1.0).
    ///   Gene pairs with relative length difference > (1 - len_dif_percent) are excluded.
    ///   Matches CD-HIT's `-s` parameter.
    ///
    /// # Returns
    ///
    /// A vector of gene clusters.
    fn cluster(&self, genes: &[Gene], identity_threshold: f32, len_dif_percent: f32) -> Result<Vec<GeneCluster>>;

    /// Get the name of this clustering backend.
    fn name(&self) -> &str;

    /// Check if this backend is available on the system.
    fn is_available(&self) -> bool;
}