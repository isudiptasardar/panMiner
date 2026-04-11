//! Core traits for downstream analysis runners.

use std::path::Path;
use crate::error::Result;

/// Input types that downstream analyses may require from the output directory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DownstreamInput {
    /// final_graph.gml - the final pangenome graph
    FinalGraph,
    /// gene_presence_absence.csv - presence/absence matrix
    PresenceAbsenceCsv,
    /// combined_protein_CDS.fasta - protein sequences for all gene clusters
    ProteinFasta,
    /// combined_DNA_CDS.fasta - DNA sequences for all gene clusters
    DnaFasta,
    /// gene_data.csv - gene cluster metadata
    GeneDataCsv,
    /// User-supplied phenotypes file (for GWAS)
    PhenotypesFile,
    /// User-supplied phylogenetic tree (Newick format)
    PhylogeneticTree,
    /// AMRFinderPlus database path
    AmrDatabase,
}

impl DownstreamInput {
    /// Returns the expected filename for this input type.
    pub fn filename(&self) -> &'static str {
        match self {
            DownstreamInput::FinalGraph => "final_graph.gml",
            DownstreamInput::PresenceAbsenceCsv => "gene_presence_absence.csv",
            DownstreamInput::ProteinFasta => "combined_protein_CDS.fasta",
            DownstreamInput::DnaFasta => "combined_DNA_CDS.fasta",
            DownstreamInput::GeneDataCsv => "gene_data.csv",
            DownstreamInput::PhenotypesFile => "phenotypes.txt",
            DownstreamInput::PhylogeneticTree => "tree.nwk",
            DownstreamInput::AmrDatabase => ".amr_db",
        }
    }
}

/// Result trait for downstream analysis outputs.
///
/// All downstream analysis runners must implement this trait to allow
/// their outputs to be written to disk and summarized.
pub trait DownstreamResult: Send + Sync {
    /// Write all output files to the specified directory.
    fn write_to(&self, dir: &Path) -> Result<()>;

    /// Summary statistics as a human-readable string.
    fn summary(&self) -> String;
}

/// Trait for all downstream analysis runners.
///
/// Implementors must be thread-safe (Send + Sync) to support parallel
/// execution of multiple downstream analyses.
pub trait DownstreamRunner: Send + Sync {
    /// Output type produced by this runner.
    type Output: DownstreamResult;

    /// Run the analysis given an output directory containing PanMiner outputs.
    ///
    /// Runners read what they need from the output directory (graph GML,
    /// P/A CSV, protein FASTA, etc.) rather than requiring in-memory graph.
    fn run(&self, output_dir: &Path) -> Result<Self::Output>;

    /// Name of the tool/analysis.
    fn name(&self) -> &str;

    /// Check if the external tool is installed and available.
    fn is_available(&self) -> bool;

    /// Declare required input files that must exist in the output directory.
    fn required_inputs(&self) -> Vec<DownstreamInput>;
}
