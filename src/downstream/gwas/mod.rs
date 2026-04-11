//! GWAS runners for gene-phenotype association analysis.
//!
//! This module provides runners for multiple GWAS tools:
//! - **pyseer** - Standard gene-based GWAS with lineage structuring
//! - **Scoary2** - Gene-trait association with support for continuous phenotypes
//! - **SpydrPick** - Co-selection and epistasis detection

pub mod pyseer;
pub mod scoary;
pub mod spydrpick;

// Re-exports
pub use pyseer::{PyseerRunner, PyseerGWASResult};
pub use scoary::{Scoary2Runner, Scoary2Result, ScoaryAssociation};
pub use spydrpick::{SpydrPickRunner, SpydrPickResult, SpydrPickCorrelation};
