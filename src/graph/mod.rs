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
mod structural_variants;
mod merge;

pub use types::{Gene, GeneId, GenomeId, Sequence, Strand, Node, Edge, EdgeKey};
pub use types::{GeneCluster, ClusterId, PangenomeGraph, GenomeMetadata};
pub use concurrent::ConcurrentGraph;
pub use matrix::BitPackedMatrix;
pub use builder::GraphBuilder;
pub use structural_variants::{StructuralVariantDetector, StructuralVariant, VariantType};
pub use merge::{merge_pangenomes, MergeResult};