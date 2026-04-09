//! CPU-based gene clustering fallback.

use rayon::prelude::*;

use crate::error::Result;
use crate::graph::{Gene, GeneCluster};
use super::traits::Clusterer;

/// CPU-based gene clusterer using greedy incremental clustering.
///
/// This is a fallback when MMseqs2 is not available. It uses a simple
/// greedy approach with Rayon for parallelism.
pub struct CpuClusterer {
    /// Number of threads to use
    threads: usize,
}

impl CpuClusterer {
    /// Create a new CPU clusterer.
    pub fn new(threads: usize) -> Self {
        Self { threads }
    }

    /// Compute sequence identity between two sequences.
pub fn simd_sequence_identity(a: &[u8], b: &[u8]) -> f32 {
    let min_len = a.len().min(b.len());
    if min_len == 0 {
        return 0.0;
    }

    let matches: usize = a[..min_len]
        .iter()
        .zip(&b[..min_len])
        .map(|(x, y)| if x == y { 1 } else { 0 })
        .sum();

    matches as f32 / min_len as f32
}


    /// Greedy incremental clustering.
    ///
    /// For each gene, compare against existing cluster centroids.
    /// If identity >= threshold, add to cluster. Otherwise, start new cluster.
    fn greedy_cluster(
        &self,
        genes: &[Gene],
        identity_threshold: f32,
    ) -> Vec<GeneCluster> {
        let mut clusters: Vec<GeneCluster> = Vec::new();
        let mut centroids: Vec<&[u8]> = Vec::new();

        for gene in genes {
            if gene.sequence.is_empty() {
                continue;
            }

            // Find best matching centroid
            let best_match = centroids
                .iter()
                .enumerate()
                .map(|(i, centroid)| (i, Self::simd_sequence_identity(centroid, &gene.sequence)))
                .filter(|(_, identity)| *identity >= identity_threshold)
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

            match best_match {
                Some((cluster_idx, _)) => {
                    clusters[cluster_idx].add_gene(gene.id.clone());
                    clusters[cluster_idx].support += 1;
                }
                None => {
                    let mut new_cluster = GeneCluster::new(format!("cluster_{}", clusters.len()));
                    new_cluster.centroid = Some(gene.sequence.clone());
                    new_cluster.add_gene(gene.id.clone());
                    new_cluster.support = 1;
                    centroids.push(&gene.sequence);
                    clusters.push(new_cluster);
                }
            }
        }

        clusters
    }
}

impl Default for CpuClusterer {
    fn default() -> Self {
        Self::new(
            std::thread::available_parallelism()
                .map(|p| p.get())
                .unwrap_or(1)
        )
    }
}

impl Clusterer for CpuClusterer {
    fn cluster(&self, genes: &[Gene], identity_threshold: f32) -> Result<Vec<GeneCluster>> {
        // Configure rayon thread pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .map_err(|e| crate::error::Error::Parallel(e.to_string()))?;

        let clusters = pool.install(|| {
            self.greedy_cluster(genes, identity_threshold)
        });

        Ok(clusters)
    }

    fn name(&self) -> &str {
        "CPU-Greedy"
    }

    fn is_available(&self) -> bool {
        true // Always available
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{GenomeId, Strand, GeneId};

    fn make_gene(id: &str, seq: &[u8]) -> Gene {
        Gene {
            id: GeneId::new(id),
            sequence: seq.to_vec(),
            genome_id: GenomeId::new("test"),
            contig: "contig1".to_string(),
            start: 0,
            end: seq.len(),
            strand: Strand::Plus,
            annotation: None,
        }
    }


    #[test]
    fn test_cpu_clustering() {
        let clusterer = CpuClusterer::new(1);
        let genes = vec![
            make_gene("g1", b"ATCGATCGATCG"),
            make_gene("g2", b"ATCGATCGATCG"), // Identical to g1
            make_gene("g3", b"NNNNNNNNNNNN"), // Very different
        ];

        let clusters = clusterer.cluster(&genes, 0.98).unwrap();

        // g1 and g2 should cluster together, g3 separate
        assert_eq!(clusters.len(), 2);
    }
}

#[cfg(test)]
mod simd_tests {
    use super::*;

    #[test]
    fn test_simd_sequence_identity() {
        let seq_a = b"ATCGATCGATCGATCGATCGATCGATCGATCG";
        let seq_b = b"ATCGATCGATCGATCGATCGATCGATCGATCG";
        let seq_c = b"ATCGATCGATCGATCGATCGATCGATCGAACG"; // 1 mismatch

        assert_eq!(CpuClusterer::simd_sequence_identity(seq_a, seq_b), 1.0);
        assert_eq!(CpuClusterer::simd_sequence_identity(seq_a, seq_c), 31.0 / 32.0);

        // Empty
        assert_eq!(CpuClusterer::simd_sequence_identity(b"", b"A"), 0.0);

        // Different lengths
        assert_eq!(CpuClusterer::simd_sequence_identity(b"ATCG", b"ATC"), 1.0); // limited by min_len
    }
}


