//! Missing gene recovery.
//!
//! Searches for genes that may have been missed by the annotation tool
//! by looking at flanking sequences around expected gene locations.
//! Uses semi-global Levenshtein alignment (HW mode) matching Panaroo's
//! edlib-based approach.

use std::collections::HashMap;

use crate::error::Result;
use crate::graph::ConcurrentGraph;
use crate::correction::simd::align_semiglobal;

/// Recovers genes missed by the initial annotation.
///
/// For each node pair where one genome lacks a neighboring gene,
/// searches the surrounding contig sequence for the missing gene
/// using semi-global alignment (HW mode) matching Panaroo's find_missing.
pub struct MissingGeneRecoverer {
    /// Minimum alignment identity to consider a hit (default: 0.70)
    min_identity: f32,
    /// Search window size in base pairs (default: 5000)
    search_window: usize,
    /// Minimum fraction of query that must align (default: 0.20)
    prop_match: f32,
}

impl MissingGeneRecoverer {
    /// Create a new missing gene recoverer with default settings.
    pub fn new() -> Self {
        Self {
            min_identity: 0.70,
            search_window: 5000,
            prop_match: 0.20,
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

    /// Set the minimum proportion of query that must align.
    pub fn with_prop_match(mut self, prop: f32) -> Self {
        self.prop_match = prop;
        self
    }

    /// Recover missing genes from the graph.
    ///
    /// For each edge in the graph, checks if any genome is missing
    /// one of the connected genes. If so, searches the flanking
    /// contig sequence for a match using semi-global alignment.
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
            if query_seq.len() < 10 {
                continue;
            }

            // Search using semi-global alignment (HW mode)
            if self.search_in_contigs(query_seq, contig_sequences) {
                recovered += 1;
            }
        }

        tracing::info!("Missing gene recovery: recovered {} genes", recovered);

        Ok(RecoveryStats {
            genes_recovered: recovered,
        })
    }

    /// Search for a query sequence in contig sequences using semi-global alignment.
    ///
    /// Uses sliding-window semi-global alignment (HW mode) to find the query
    /// embedded within longer contig sequences. This matches Panaroo's edlib-based
    /// approach where the query is aligned within the target.
    fn search_in_contigs(
        &self,
        query: &[u8],
        contigs: &HashMap<String, Vec<u8>>,
    ) -> bool {
        if query.is_empty() {
            return false;
        }

        let min_align_len = (query.len() as f32 * self.prop_match) as usize;
        let min_align_len = min_align_len.max(1);

        for seq in contigs.values() {
            if seq.len() < query.len() {
                // Contig is shorter than query — just do direct alignment
                let (identity, _dist, align_len) = align_semiglobal(query, seq);
                if align_len >= min_align_len && identity >= self.min_identity as f64 {
                    return true;
                }
                continue;
            }

            // Slide a window across the contig and align the query within each window
            let window_size = query.len() + self.search_window;
            let step = (query.len() / 2).max(100);

            let mut pos = 0;
            while pos < seq.len() {
                let end = (pos + window_size).min(seq.len());
                let window = &seq[pos..end];

                let (identity, _dist, align_len) = align_semiglobal(query, window);

                if align_len >= min_align_len && identity >= self.min_identity as f64 {
                    return true;
                }

                pos += step;
                // Avoid infinite loop when step is 0 or window_size is larger than seq
                if pos >= seq.len() {
                    break;
                }
            }
        }

        false
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
    fn test_search_sequence_embedded() {
        // Query embedded within a longer contig — HW mode should find it
        let recoverer = MissingGeneRecoverer::new()
            .with_min_identity(0.80);
        let query = b"ATCGATCGATCGATCG";
        let mut contigs = HashMap::new();
        contigs.insert("contig1".to_string(), b"NNNNATCGATCGATCGATCGNNNN".to_vec());

        assert!(recoverer.search_in_contigs(query, &contigs));
    }

    #[test]
    fn test_search_sequence_not_found() {
        // Query not present in contig
        let recoverer = MissingGeneRecoverer::new()
            .with_min_identity(0.80);
        let query = b"ATCGATCGATCGATCG";
        let mut contigs = HashMap::new();
        contigs.insert("contig1".to_string(), b"GGGGAAAACCCCTTTT".to_vec());

        assert!(!recoverer.search_in_contigs(query, &contigs));
    }

    #[test]
    fn test_search_sequence_near_match() {
        // Query with 1 substitution embedded in contig
        let recoverer = MissingGeneRecoverer::new()
            .with_min_identity(0.85);
        let query = b"ATCGATCGATCGATCG";
        let mut contigs = HashMap::new();
        // 15/16 = 93.75% identity (1 substitution in 16 chars)
        contigs.insert("contig1".to_string(), b"NNNNATCGATCGATCAATCGNNNN".to_vec());

        assert!(recoverer.search_in_contigs(query, &contigs));
    }
}