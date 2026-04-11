//! PanMiner - A modern pangenome analysis tool with GPU and CPU support.
//!
//! This library provides efficient pangenome construction and analysis
//! capabilities for prokaryotic genomes, implementing a graph-based approach
//! similar to Panaroo with modern optimizations.
//!
//! # Features
//!
//! - Memory-mapped I/O for large datasets
//! - Parallel processing with Rayon
//! - GPU-accelerated clustering via MMseqs2
//! - Concurrent graph construction with DashMap
//! - Multiple output formats (CSV, FASTA, GML, JSON, Parquet, HTML)
//!
//! # Example
//!
//! ```no_run
//! use panminer::{PanminerConfig, PanminerPipeline};
//!
//! # fn main() -> panminer::Result<()> {
//! let config = PanminerConfig::default();
//! let pipeline = PanminerPipeline::new(config);
//! let result = pipeline.run()?;
//! # Ok(())
//! # }
//! ```

#[allow(dead_code)]
pub mod config;
pub mod error;
pub mod io;
pub mod clustering;
pub mod graph;
pub mod correction;
pub mod output;
pub mod gwas;
pub mod downstream;
pub mod pipeline;

// Re-exports for convenience
pub use config::PanminerConfig;
pub use error::{Error, Result};
pub use pipeline::PanminerPipeline;
pub use graph::{PangenomeGraph, Gene, GeneCluster, Node, Edge, merge_pangenomes, MergeResult};
pub use config::OutputFormat;
pub use io::{BaktaRunner, BaktaDbType, translate, translate_with_stop};

/// Prelude module for common imports
pub mod prelude {
    pub use crate::config::PanminerConfig;
    pub use crate::error::{Error, Result};
    pub use crate::pipeline::PanminerPipeline;
    pub use crate::graph::{PangenomeGraph, Gene, GeneCluster};
}

/// Library version
pub const VERSION: &str = env!("CARGO_PKG_VERSION");