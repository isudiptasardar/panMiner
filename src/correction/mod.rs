//! Error correction module for PanMiner.
//!
//! Implements graph-based error correction steps:
//! - Contamination removal
//! - Fragment merging (mistranslation correction)
//! - Missing gene recovery

mod contamination;
mod fragment;
mod missing;

pub use contamination::ContaminationRemover;
pub use fragment::{FragmentMerger, MergeStats};
pub use missing::{MissingGeneRecoverer, RecoveryStats};