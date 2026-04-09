//! Bit-packed presence/absence matrix for memory efficiency.

use std::path::Path;
use csv::Writer;

use crate::error::Result;

/// Bit-packed presence/absence matrix.
///
/// This matrix stores gene presence/absence in a memory-efficient format
/// using 1 bit per genome per cluster, providing 8x memory reduction
/// compared to a naive boolean matrix.
///
/// # Memory Usage
///
/// For 1000 genomes and 10,000 clusters:
/// - Naive storage: 1000 * 10000 = 10,000,000 bytes
/// - Bit-packed: 10,000,000 / 8 = 1,250,000 bytes
pub struct BitPackedMatrix {
    /// Packed data (8 genomes per byte)
    data: Vec<u8>,
    /// Number of genomes
    num_genomes: usize,
    /// Number of clusters
    num_clusters: usize,
    /// Bytes per row (num_genomes rounded up to nearest 8)
    bytes_per_row: usize,
    /// Genome names (column headers)
    genome_names: Vec<String>,
    /// Cluster IDs (row identifiers)
    cluster_ids: Vec<String>,
}

impl BitPackedMatrix {
    /// Create a new bit-packed matrix.
    ///
    /// # Arguments
    ///
    /// * `num_genomes` - Number of genomes (columns)
    /// * `num_clusters` - Number of clusters (rows)
    pub fn new(num_genomes: usize, num_clusters: usize) -> Self {
        let bytes_per_row = (num_genomes + 7) / 8;
        let data = vec![0u8; num_clusters * bytes_per_row];

        Self {
            data,
            num_genomes,
            num_clusters,
            bytes_per_row,
            genome_names: Vec::with_capacity(num_genomes),
            cluster_ids: Vec::with_capacity(num_clusters),
        }
    }

    /// Set genome names (column headers).
    pub fn set_genome_names(&mut self, names: Vec<String>) {
        self.genome_names = names;
    }

    /// Set cluster IDs (row identifiers).
    pub fn set_cluster_ids(&mut self, ids: Vec<String>) {
        self.cluster_ids = ids;
    }

    /// Set the presence/absence value for a cell.
    ///
    /// # Arguments
    ///
    /// * `genome` - Genome index (column)
    /// * `cluster` - Cluster index (row)
    /// * `present` - Whether the gene is present
    #[inline]
    pub fn set(&mut self, genome: usize, cluster: usize, present: bool) {
        debug_assert!(genome < self.num_genomes);
        debug_assert!(cluster < self.num_clusters);

        let byte_idx = genome / 8;
        let bit_idx = genome % 8;
        let row_start = cluster * self.bytes_per_row;

        if present {
            self.data[row_start + byte_idx] |= 1 << bit_idx;
        } else {
            self.data[row_start + byte_idx] &= !(1 << bit_idx);
        }
    }

    /// Get the presence/absence value for a cell.
    ///
    /// # Arguments
    ///
    /// * `genome` - Genome index (column)
    /// * `cluster` - Cluster index (row)
    #[inline]
    pub fn get(&self, genome: usize, cluster: usize) -> bool {
        debug_assert!(genome < self.num_genomes);
        debug_assert!(cluster < self.num_clusters);

        let byte_idx = genome / 8;
        let bit_idx = genome % 8;
        let row_start = cluster * self.bytes_per_row;

        (self.data[row_start + byte_idx] >> bit_idx) & 1 == 1
    }

    /// Count how many genomes have this cluster (row sum).
    pub fn count_present(&self, cluster: usize) -> usize {
        debug_assert!(cluster < self.num_clusters);

        let row_start = cluster * self.bytes_per_row;
        let row = &self.data[row_start..row_start + self.bytes_per_row];

        row.iter().map(|&b| b.count_ones() as usize).sum()
    }

    /// Count how many clusters a genome has (column sum).
    pub fn count_genome_clusters(&self, genome: usize) -> usize {
        debug_assert!(genome < self.num_genomes);

        let bit_idx = genome % 8;
        let mask = 1u8 << bit_idx;

        (0..self.num_clusters)
            .filter(|&cluster| {
                let row_start = cluster * self.bytes_per_row;
                let byte_idx = genome / 8;
                (self.data[row_start + byte_idx] & mask) != 0
            })
            .count()
    }

    /// Get the total number of genes present in the matrix.
    pub fn total_present(&self) -> usize {
        self.data.iter().map(|&b| b.count_ones() as usize).sum()
    }

    /// Check if a cluster is core (present in all genomes).
    pub fn is_core(&self, cluster: usize, threshold: f32) -> bool {
        let count = self.count_present(cluster);
        (count as f32 / self.num_genomes as f32) >= threshold
    }

    /// Check if a cluster is accessory.
    pub fn is_accessory(&self, cluster: usize, core_threshold: f32) -> bool {
        !self.is_core(cluster, core_threshold)
    }

    /// Export to CSV format.
    pub fn to_csv(&self, path: &Path) -> Result<()> {
        let mut writer = Writer::from_path(path)?;

        // Header row
        let mut header = vec!["Gene".to_string()];
        header.extend(self.genome_names.iter().cloned());
        writer.write_record(&header)?;

        // Data rows
        for (cluster_idx, cluster_id) in self.cluster_ids.iter().enumerate() {
            let mut row = vec![cluster_id.clone()];
            for genome_idx in 0..self.num_genomes {
                row.push(if self.get(genome_idx, cluster_idx) {
                    "1".to_string()
                } else {
                    "0".to_string()
                });
            }
            writer.write_record(&row)?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Export to TSV format (tab-separated).
    pub fn to_tsv(&self, path: &Path) -> Result<()> {
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(path)?;

        // Header row
        let mut header = vec!["Gene".to_string()];
        header.extend(self.genome_names.iter().cloned());
        writer.write_record(&header)?;

        // Data rows
        for (cluster_idx, cluster_id) in self.cluster_ids.iter().enumerate() {
            let mut row = vec![cluster_id.clone()];
            for genome_idx in 0..self.num_genomes {
                row.push(if self.get(genome_idx, cluster_idx) {
                    "1".to_string()
                } else {
                    "0".to_string()
                });
            }
            writer.write_record(&row)?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Get the memory usage in bytes.
    pub fn memory_usage(&self) -> usize {
        self.data.len()
    }

    /// Get the number of genomes.
    pub fn num_genomes(&self) -> usize {
        self.num_genomes
    }

    /// Get the number of clusters.
    pub fn num_clusters(&self) -> usize {
        self.num_clusters
    }
}

use csv::WriterBuilder;

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_bit_packed_matrix_basic() {
        let mut matrix = BitPackedMatrix::new(10, 5);

        // Set some values
        matrix.set(0, 0, true);
        matrix.set(1, 0, true);
        matrix.set(0, 1, false);

        // Check values
        assert!(matrix.get(0, 0));
        assert!(matrix.get(1, 0));
        assert!(!matrix.get(0, 1));

        // Count present
        assert_eq!(matrix.count_present(0), 2);
    }

    #[test]
    fn test_bit_packed_matrix_memory() {
        let matrix = BitPackedMatrix::new(1000, 10000);

        // Should use about 10000 * 125 bytes = 1.25 MB
        let usage = matrix.memory_usage();
        assert!(usage < 2_000_000); // Less than 2 MB
    }

    #[test]
    fn test_bit_packed_matrix_core_detection() {
        let mut matrix = BitPackedMatrix::new(10, 3);

        // Cluster 0: present in all genomes (core)
        for g in 0..10 {
            matrix.set(g, 0, true);
        }

        // Cluster 1: present in 5 genomes (accessory)
        for g in 0..5 {
            matrix.set(g, 1, true);
        }

        assert!(matrix.is_core(0, 0.95));
        assert!(!matrix.is_core(1, 0.95));
        assert_eq!(matrix.count_present(0), 10);
        assert_eq!(matrix.count_present(1), 5);
    }

    #[test]
    fn test_bit_packed_matrix_csv_export() {
        let mut matrix = BitPackedMatrix::new(3, 2);
        matrix.set_genome_names(vec!["g1".to_string(), "g2".to_string(), "g3".to_string()]);
        matrix.set_cluster_ids(vec!["c1".to_string(), "c2".to_string()]);

        matrix.set(0, 0, true);
        matrix.set(1, 0, true);
        matrix.set(2, 1, true);

        let temp = NamedTempFile::new().unwrap();
        matrix.to_csv(temp.path()).unwrap();

        let content = std::fs::read_to_string(temp.path()).unwrap();
        assert!(content.contains("Gene"));
        assert!(content.contains("c1"));
        assert!(content.contains("g1"));
    }
}