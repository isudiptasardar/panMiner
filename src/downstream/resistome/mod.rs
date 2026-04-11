//! Resistome/AMR detection runners.
//!
//! This module provides runners for antimicrobial resistance (AMR) gene
//! detection using AMRFinderPlus, which offers NCBI-curated detection with
//! hierarchical evidence levels (EXACT > ALLELE > BLAST > HMM > PARTIAL).
//!
//! # Features
//!
//! - Hierarchical detection method with evidence tracking
//! - Taxon-specific point mutation analysis
//! - Stress response and virulence factor detection
//! - Continuous database updates via NCBI
//!
//! # Reference
//!
//! Feldgarden et al. (2021) "AMRFinderPlus: a curated and benchmarked
//! tool for detection of AMR genes in bacterial nucleotide sequences."
//! *Scientific Reports* 11:12720.
//!
//! # Example
//!
//! ```no_run
//! use std::path::Path;
//! use panminer::downstream::{DownstreamRunner, DownstreamResult, AmrFinderRunner};
//!
//! # fn main() -> panminer::Result<()> {
//! if let Some(runner) = AmrFinderRunner::detect() {
//!     let output = runner.run(Path::new("panminer_output"))?;
//!     output.write_to(Path::new("downstream_results"))?;
//! }
//! # Ok(())
//! # }
//! ```

//! Resistome/AMR detection runners.
pub mod amrfinder;

// Re-exports
pub use amrfinder::{AmrFinderRunner, AmrFinderResult, AmrGene};
