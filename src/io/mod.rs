//! I/O module for PanMiner.
//!
//! This module provides memory-mapped file parsing for GFF3 and FASTA formats,
//! compressed intermediate storage, and streaming pipeline support.

mod mmap;
mod gff;
mod fasta;
mod compress;
mod streaming;
mod qc_traits;

pub use mmap::MmapFile;
pub use gff::GffParser;
pub use fasta::FastaParser;
pub use compress::{CompressedWriter, CompressedReader};
pub use streaming::{StreamingPipeline, PartialGraph};
pub use qc_traits::{GenomeQC, QcMode, QcRunner, CheckmQcRunner};