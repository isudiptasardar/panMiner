//! Evolutionary model runners for gene gain/loss rate estimation.
//!
//! This module provides runners for evolutionary analysis tools:
//! - **Panstripe** - Phylogenetically-informed GLM for gene gain/loss rates
//!
//! Panstripe estimates gene gain and loss rates using generalized linear
//! models, providing superior robustness to annotation errors compared to
//! traditional IMG/FMG approaches.
//!
//! # Reference
//!
//! Tonkin-Hill et al. (2023) "Panstripe: phylogenetically-informed gene
//! gain and loss rates with lineage-specific detection." *Genome Research*.
//!
//! # Example
//!
//! ```no_run
//! use std::path::Path;
//! use panminer::downstream::{DownstreamRunner, DownstreamResult, PanstripeRunner};
//!
//! # fn main() -> panminer::Result<()> {
//! let runner = PanstripeRunner::new();
//! if runner.is_available() {
//!     let output = runner.run(Path::new("panminer_output"))?;
//!     output.write_to(Path::new("downstream_results"))?;
//! }
//! # Ok(())
//! # }
//! ```

//! Evolutionary model runners.
pub mod panstripe;
pub mod pangrowth;

// Re-exports
pub use panstripe::{PanstripeRunner, PanstripeResult};
pub use pangrowth::{PangrowthRunner, PangrowthResult, OpennessClassification};
