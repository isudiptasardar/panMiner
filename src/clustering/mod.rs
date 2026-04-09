//! Clustering module for PanMiner.
//!
//! This module provides gene clustering functionality using MMseqs2-GPU
//! or a CPU fallback.

mod traits;
mod mmseqs;
mod cpu;
mod alignment_traits;
mod alignment_mafft;
mod alignment_clustal;
mod alignment_prank;

pub use traits::Clusterer;
pub use mmseqs::MMseqsRunner;
pub use cpu::CpuClusterer;
pub use alignment_traits::{AlignmentRunner, AlignmentTool, AlignmentResult, build_alignment_from_graph};
pub use alignment_mafft::MafftRunner;
pub use alignment_clustal::ClustalOmegaRunner;
pub use alignment_prank::PrankRunner;