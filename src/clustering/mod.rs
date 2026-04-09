//! Clustering module for PanMiner.
//!
//! This module provides gene clustering functionality using MMseqs2-GPU
//! or a CPU fallback.

mod traits;
mod mmseqs;
mod cpu;

pub use traits::Clusterer;
pub use mmseqs::MMseqsRunner;
pub use cpu::CpuClusterer;