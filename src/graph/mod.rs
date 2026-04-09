//! Graph module for PanMiner.
//!
//! This module provides the core data structures for the pangenome graph:
//! - Gene: A gene from a genome
//! - GeneCluster: A cluster of orthologous genes
//! - PangenomeGraph: The main graph structure
//!
//! The graph uses DashMap for concurrent updates and supports parallel processing.

mod types;
mod concurrent;
mod matrix;
mod builder;

pub use types::{Gene, GeneId, GenomeId, Sequence, Strand, Node, Edge, EdgeKey};
pub use types::{GeneCluster, ClusterId, PangenomeGraph};
pub use concurrent::ConcurrentGraph;
pub use matrix::BitPackedMatrix;
pub use builder::GraphBuilder;