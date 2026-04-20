//! Downstream analysis module for PanMiner.
//!
//! This module provides downstream analysis capabilities including GWAS,
//! evolutionary modeling, AMR detection, and pangenome exploration tools.
//!
//! # Architecture
//!
//! Each analysis runner implements the [`DownstreamRunner`] trait, which
//! operates on PanMiner output directories rather than requiring in-memory
//! graph access. This allows users to run analyses independently after
//! pangenome construction.
//!
//! # Module Structure
//!
//! - [`gwas`] - GWAS runners (pyseer, Scoary2, SpydrPick)
//! - [`evolution`] - Evolutionary model runners (Panstripe)
//! - [`resistome`] - AMR detection runners (AMRFinderPlus)
//! - [`exploration`] - Exploration and visualization tools (neighborhood extraction, accumulation curves)

pub mod traits;
pub mod gwas;
pub mod evolution;
pub mod resistome;
pub mod exploration;

// Re-exports
pub use traits::{DownstreamInput, DownstreamResult, DownstreamRunner};

// Re-export all runners
pub use gwas::{PyseerRunner, Scoary2Runner, SpydrPickRunner};
pub use evolution::{PanstripeRunner, PangrowthRunner};
pub use resistome::AmrFinderRunner;
pub use exploration::{GeneNeighborhoodExtractor, AccumulationCurveRunner, GrapeTreeExportRunner, ItolAnnotationRunner};
