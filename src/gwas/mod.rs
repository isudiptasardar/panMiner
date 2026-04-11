//! GWAS (Genome-Wide Association Studies) module.
//!
//! Provides integration with GWAS tools (pyseer) for association
//! analysis on pangenome data.

pub mod traits;
pub mod pyseer;

pub use traits::{GWASRunner, GWASOutput, GWASResult};
pub use pyseer::PyseerRunner;