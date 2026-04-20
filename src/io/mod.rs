//! I/O module for PanMiner.
//!
//! This module provides memory-mapped file parsing for GFF3 and FASTA formats,
//! compressed intermediate storage, streaming pipeline support, and QC tool
//! integration (CheckM2, skani).

mod mmap;
mod gff;
mod fasta;
mod compress;
mod streaming;
mod qc_traits;
mod bakta;
mod translate;
mod skani;
mod mds;
mod ggcat;
mod ggcaller;
mod subprocess;
mod integrate;
pub mod extract_gene;

#[cfg(feature = "prodigal")]
pub mod orphos;

pub use mmap::MmapFile;
pub use extract_gene::extract_gene;
pub use gff::GffParser;
pub use fasta::FastaParser;
pub use compress::{CompressedWriter, CompressedReader};
pub use streaming::{StreamingPipeline, PartialGraph};
pub use qc_traits::{GenomeQC, QcMode, QcRunner, CheckmQcRunner};
pub use bakta::{BaktaRunner, BaktaDbType, is_gff_file, is_genbank_file, genbank_to_fasta};
pub use translate::{translate, translate_with_stop};
pub use skani::SkaniRunner;
pub use mds::{compute_mds, compute_mds_with_labels, MdsProjection};
pub use ggcat::{GGCATBuilder, CDBGGraph, CDBGStats, compute_cdbg_stats};
pub use ggcaller::{GGCallerRunner, GGCallerOutput};
pub use subprocess::run_with_timeout;
pub use integrate::{integrate_genome, IntegrateResult};

#[cfg(feature = "prodigal")]
pub use orphos::{OrphosRunner, PredictedGene, Strand};