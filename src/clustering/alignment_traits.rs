//! Alignment runner trait and implementations.
//!
//! Provides trait for running multiple sequence alignment (MSA) tools.
//! Supports MAFFT, Clustal Omega, and PRANK.

use crate::error::Result;

/// Result of running an MSA tool.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Number of sequences aligned
    pub num_sequences: usize,
    /// Length of alignment (columns)
    pub alignment_length: usize,
    /// Aligned sequences in FASTA format
    pub aligned_fasta: String,
    /// Tool-specific metadata
    pub tool: AlignmentTool,
}

/// Multiple sequence alignment tools.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentTool {
    /// MAFFT - Fast multiple sequence alignment
    Mafft,
    /// Clustal Omega - Progressive multiple sequence alignment
    ClustalOmega,
    /// PRANK - Phylogeny-aware multiple sequence alignment
    Prank,
}

impl AlignmentTool {
    /// Get the tool name for CLI invocation.
    pub fn name(&self) -> &'static str {
        match self {
            AlignmentTool::Mafft => "MAFFT",
            AlignmentTool::ClustalOmega => "Clustal Omega",
            AlignmentTool::Prank => "PRANK",
        }
    }

    /// Get the executable name for this tool.
    pub fn executable(&self) -> &'static str {
        match self {
            AlignmentTool::Mafft => "mafft",
            AlignmentTool::ClustalOmega => "clustalo",
            AlignmentTool::Prank => "prank",
        }
    }
}

/// Trait for MSA runner implementations.
pub trait AlignmentRunner {
    /// Run multiple sequence alignment on the given sequences.
    ///
    /// # Arguments
    ///
    /// * `sequences` - Vector of (name, sequence) pairs
    /// * `tool` - The alignment tool to use
    ///
    /// # Returns
    ///
    /// AlignmentResult containing aligned sequences
    fn run_msa(&self, sequences: &[(String, Vec<u8>)], tool: AlignmentTool) -> Result<AlignmentResult>;

    /// Get the name of this alignment runner.
    fn name(&self) -> &str;

    /// Check if this runner is available on the system.
    fn is_available(&self) -> bool;
}

/// Build an alignment from a PangenomeGraph.
pub fn build_alignment_from_graph(
    graph: &crate::graph::PangenomeGraph,
) -> Vec<(String, Vec<u8>)> {
    // Collect sequences from nodes that have centroids
    let mut sequences = Vec::new();

    for (cluster_id, node) in &graph.nodes {
        if let Some(centroid) = node.centroid_sequences.first() {
            sequences.push((cluster_id.to_string(), centroid.clone()));
        }
    }

    sequences
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alignment_tool_names() {
        assert_eq!(AlignmentTool::Mafft.name(), "MAFFT");
        assert_eq!(AlignmentTool::ClustalOmega.name(), "Clustal Omega");
        assert_eq!(AlignmentTool::Prank.name(), "PRANK");
    }

    #[test]
    fn test_alignment_tool_executables() {
        assert_eq!(AlignmentTool::Mafft.executable(), "mafft");
        assert_eq!(AlignmentTool::ClustalOmega.executable(), "clustalo");
        assert_eq!(AlignmentTool::Prank.executable(), "prank");
    }
}
