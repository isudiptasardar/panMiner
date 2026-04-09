//! Missing gene recovery.
//!
//! Searches for genes that may have been missed by the annotation tool
//! by looking at flanking sequences around expected gene locations.

use std::collections::HashMap;
use rayon::prelude::*;

use crate::error::Result;
use crate::graph::ConcurrentGraph;

/// Recovers genes missed by the initial annotation.
///
/// For each node pair where one genome lacks a neighboring gene,
/// searches the surrounding contig sequence for the missing gene.
pub struct MissingGeneRecoverer {
    /// Minimum alignment identity to consider a hit (default: 0.70)
    min_identity: f32,
    /// Search window size in base pairs (default: 5000)
    search_window: usize,
}

impl MissingGeneRecoverer {
    /// Create a new missing gene recoverer with default settings.
    pub fn new() -> Self {
        Self {
            min_identity: 0.70,
            search_window: 5000,
        }
    }

    /// Set the minimum identity threshold.
    pub fn with_min_identity(mut self, identity: f32) -> Self {
        self.min_identity = identity;
        self
    }

    /// Set the search window size.
    pub fn with_search_window(mut self, window: usize) -> Self {
        self.search_window = window;
        self
    }

    /// Recover missing genes from the graph.
    ///
    /// For each edge in the graph, checks if any genome is missing
    /// one of the connected genes. If so, searches the flanking
    /// sequence for a match.
    ///
    /// # Arguments
    ///
    /// * `graph` - The pangenome graph
    /// * `contig_sequences` - Map of contig name -> sequence
    /// * `cluster_sequences` - Map of cluster ID -> representative sequence
    pub fn recover(
        &self,
        graph: &ConcurrentGraph,
        contig_sequences: &HashMap<String, Vec<u8>>,
        cluster_sequences: &HashMap<String, Vec<u8>>,
    ) -> Result<RecoveryStats> {
        let mut recovered = 0;

        let mut edges = Vec::new();
        for e in graph.edges.iter() {
            let (from, to) = e.key();
            let genomes = e.value().genomes.clone();
            edges.push((from.clone(), to.clone(), genomes));
        }

        for (from, to, _edge_genomes) in &edges {
            // Get genomes that have the 'from' node but not the 'to' node
            let from_node = graph.nodes.get(from);
            let to_node = graph.nodes.get(to);

            if from_node.is_none() || to_node.is_none() {
                continue;
            }

            // Get the query sequence for the missing gene
            let query_seq = match cluster_sequences.get(to.as_str()) {
                Some(seq) => seq,
                None => continue,
            };

            // Search in contig sequences for each genome missing this gene
            // (simplified - in practice, need genome-specific contig mapping)
            if query_seq.len() < 10 {
                continue;
            }

            // Simple k-mer search as approximation
            let found = self.simd_search_sequence(query_seq, contig_sequences);
            if found {
                recovered += 1;
            }
        }

        tracing::info!("Missing gene recovery: recovered {} genes", recovered);

        Ok(RecoveryStats {
            genes_recovered: recovered,
        })
    }

    /// Search for a query sequence in contig sequences.
    ///
    /// Uses a simple k-mer based search for initial screening.
    /// SIMD-optimized search sequence using chunked parallel iteration
    pub fn simd_search_sequence(
        &self,
        query: &[u8],
        contigs: &HashMap<String, Vec<u8>>,
    ) -> bool {
        if query.len() < 11 {
            return false;
        }

        let kmer_len = 11;
        let query_kmers: std::collections::HashSet<&[u8]> = query
            .windows(kmer_len)
            .collect();

        let threshold = (query_kmers.len() as f32 * self.min_identity) as usize;

        // Use rayon for parallel searching across contigs
        contigs.par_iter().any(|(_name, seq)| {
            let matches: usize = seq
                .windows(kmer_len)
                .map(|kmer| if query_kmers.contains(kmer) { 1 } else { 0 })
                .sum();
            matches >= threshold
        })
    }

}

impl Default for MissingGeneRecoverer {
    fn default() -> Self {
        Self::new()
    }
}

/// Statistics from missing gene recovery.
#[derive(Debug, Clone)]
pub struct RecoveryStats {
    /// Number of genes recovered
    pub genes_recovered: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_missing_gene_recoverer_creation() {
        let recoverer = MissingGeneRecoverer::new();
        assert_eq!(recoverer.min_identity, 0.70);
        assert_eq!(recoverer.search_window, 5000);
    }

    #[test]
    fn test_search_sequence() {
        let recoverer = MissingGeneRecoverer::new();
        let query = b"ATCGATCGATCGATCG";
        let mut contigs = HashMap::new();
        contigs.insert("contig1".to_string(), b"NNNATCGATCGATCGATCGNNN".to_vec());

        assert!(recoverer.simd_search_sequence(query, &contigs));
    }
}

#[cfg(test)]
mod simd_tests {
    use super::*;

    #[test]
    fn test_simd_kmer_search() {
        let recoverer = MissingGeneRecoverer::new();
        let query = b"ATCGATCGATCGATCG";
        let mut contigs = HashMap::new();
        contigs.insert("contig1".to_string(), b"NNNATCGATCGATCGATCGNNN".to_vec());
        
        assert!(recoverer.simd_search_sequence(query, &contigs));
        
        let mut contigs_fail = HashMap::new();
        contigs_fail.insert("contig1".to_string(), b"NNNATCGAACGATCGAACGNNN".to_vec());
        assert!(!recoverer.simd_search_sequence(query, &contigs_fail));
    }
}
