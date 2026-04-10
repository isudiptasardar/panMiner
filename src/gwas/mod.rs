//! GWAS (Genome-Wide Association Studies) module.
//!
//! Provides integration with pyseer for GWAS analysis on pangenome data.

pub mod traits;
pub mod pyseer;

pub use traits::{GWASRunner, GWASOutput, GWASResult};
pub use pyseer::PyseerRunner;
