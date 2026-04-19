//! CPU-based gene clustering fallback.

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
    ///
    /// Uses naive byte-by-byte comparison. For large-scale clustering,
    /// MMseqs2-GPU or SIMD-accelerated paths should be preferred.
    pub fn sequence_identity(a: &[u8], b: &[u8]) -> f32 {
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
    /// Gene pairs with relative length difference > (1 - len_dif_percent) are skipped.
    fn greedy_cluster(
        &self,
        genes: &[Gene],
        identity_threshold: f32,
        len_dif_percent: f32,
    ) -> Vec<GeneCluster> {
        let mut clusters: Vec<GeneCluster> = Vec::new();
        let mut centroids: Vec<&Gene> = Vec::new();

        for gene in genes {
            if gene.sequence.is_empty() {
                continue;
            }

            // Find best matching centroid
            let best_match = centroids
                .iter()
                .enumerate()
                .filter_map(|(i, centroid)| {
                    // Length filter: skip if relative length difference exceeds (1 - len_dif_percent)
                    let max_len = centroid.sequence.len().max(gene.sequence.len()) as f32;
                    if max_len == 0.0 {
                        return None;
                    }
                    let len_diff = (centroid.sequence.len().abs_diff(gene.sequence.len())) as f32 / max_len;
                    if len_diff > (1.0 - len_dif_percent) {
                        return None;
                    }

                    let identity = Self::sequence_identity(&centroid.sequence, &gene.sequence);
                    if identity >= identity_threshold {
                        Some((i, identity))
                    } else {
                        None
                    }
                })
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

            match best_match {
                Some((cluster_idx, _)) => {
                    clusters[cluster_idx].add_gene(gene.id.clone());
                    clusters[cluster_idx].support += 1;
                }
                None => {
                    let mut new_cluster = GeneCluster::new(format!("cluster_{}", clusters.len()));
                    new_cluster.centroids = vec![gene.sequence.clone()];
                    new_cluster.add_gene(gene.id.clone());
                    new_cluster.support = 1;
                    centroids.push(gene);
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
    fn cluster(&self, genes: &[Gene], identity_threshold: f32, len_dif_percent: f32) -> Result<Vec<GeneCluster>> {
        // Configure rayon thread pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .map_err(|e| crate::error::Error::Parallel(e.to_string()))?;

        let clusters = pool.install(|| {
            self.greedy_cluster(genes, identity_threshold, len_dif_percent)
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

        let clusters = clusterer.cluster(&genes, 0.98, 0.98).unwrap();

        // g1 and g2 should cluster together, g3 separate
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn test_length_filter_rejects_different_lengths() {
        let clusterer = CpuClusterer::new(1);

        // Create a 300bp gene and a 150bp gene with matching prefixes
        let long_seq: Vec<u8> = b"ATCG".repeat(75);   // 300 bp
        let short_seq: Vec<u8> = long_seq[..150].to_vec(); // 150 bp (prefix of long_seq)

        let genes = vec![
            make_gene("g_long", &long_seq),
            make_gene("g_short", &short_seq),
        ];

        // At len_dif_percent=0.98, max allowed difference is 2%
        // 150/300 = 50% length ratio -> 50% difference -> should be rejected
        let clusters_strict = clusterer.cluster(&genes, 0.90, 0.98).unwrap();
        assert_eq!(clusters_strict.len(), 2, "Strict len_dif_percent should separate 300bp and 150bp genes");

        // At len_dif_percent=0.50, 50% difference is allowed
        let clusters_loose = clusterer.cluster(&genes, 0.90, 0.50).unwrap();
        assert_eq!(clusters_loose.len(), 1, "Loose len_dif_percent should cluster 300bp and 150bp genes together");
    }
}

#[cfg(test)]
mod sequence_identity_tests {
    use super::*;

    #[test]
    fn test_sequence_identity() {
        let seq_a = b"ATCGATCGATCGATCGATCGATCGATCGATCG";
        let seq_b = b"ATCGATCGATCGATCGATCGATCGATCGATCG";
        let seq_c = b"ATCGATCGATCGATCGATCGATCGATCGAACG"; // 1 mismatch

        assert_eq!(CpuClusterer::sequence_identity(seq_a, seq_b), 1.0);
        assert_eq!(CpuClusterer::sequence_identity(seq_a, seq_c), 31.0 / 32.0);

        // Empty
        assert_eq!(CpuClusterer::sequence_identity(b"", b"A"), 0.0);

        // Different lengths
        assert_eq!(CpuClusterer::sequence_identity(b"ATCG", b"ATC"), 1.0); // limited by min_len
    }
}


