//! Error correction module for PanMiner.
//!
//! Implements graph-based error correction steps:
//! - Paralog resolution with synteny context
//! - Contamination removal
//! - Fragment merging (mistranslation correction)
//! - Missing gene recovery
//! - Misassembly edge cleaning

mod contamination;
mod contig_end;
mod fragment;
mod missing;
mod misassembly;
mod paralog;
mod simd;

pub use contamination::ContaminationRemover;
pub use contig_end::{ContigEndPruner, PruningStats};
pub use fragment::{FragmentMerger, MergeStats, DistanceCache};
pub use missing::{MissingGeneRecoverer, RecoveryStats};
pub use misassembly::{MisassemblyEdgeCleaner, CleaningStats};
pub use paralog::{ParalogResolver, ParalogStats};
pub use simd::compare_sequences;