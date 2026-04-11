//! Summary statistics output.
//!
//! Generates gene category counts following Panaroo conventions:
//! - Core (>=99%), Soft core (95-99%), Shell (15-95%), Cloud (0-15%)

use std::path::Path;
use crate::error::Result;
use crate::graph::BitPackedMatrix;

/// Write summary statistics to a text file.
pub fn write_summary_stats(matrix: &BitPackedMatrix, path: &Path) -> Result<()> {
    let num_genomes = matrix.num_genomes();
    if num_genomes == 0 {
        return Ok(());
    }

    let mut core = 0usize;       // >= 99%
    let mut soft_core = 0usize;  // 95-99%
    let mut shell = 0usize;      // 15-95%
    let mut cloud = 0usize;      // 0-15%

    for cluster_idx in 0..matrix.num_clusters() {
        let count = matrix.count_present(cluster_idx);
        let proportion = count as f64 / num_genomes as f64;

        if proportion >= 0.99 {
            core += 1;
        } else if proportion >= 0.95 {
            soft_core += 1;
        } else if proportion >= 0.15 {
            shell += 1;
        } else {
            cloud += 1;
        }
    }

    let total = core + soft_core + shell + cloud;

    let content = format!(
        "Core genes (>= 99%): {}\nSoft core genes (95-99%): {}\nShell genes (15-95%): {}\nCloud genes (0-15%): {}\nTotal genes: {}\n",
        core, soft_core, shell, cloud, total
    );

    std::fs::write(path, content)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_summary_stats() {
        let mut matrix = BitPackedMatrix::new(10, 4);
        matrix.set_genome_names(vec![
            "g1".to_string(), "g2".to_string(), "g3".to_string(), "g4".to_string(),
            "g5".to_string(), "g6".to_string(), "g7".to_string(), "g8".to_string(),
            "g9".to_string(), "g10".to_string(),
        ]);
        matrix.set_cluster_ids(vec![
            "c1".to_string(), "c2".to_string(), "c3".to_string(), "c4".to_string(),
        ]);

        // Core (10/10 = 100%)
        for g in 0..10 { matrix.set(g, 0, true); }
        // Soft core: make it exactly 95% (9.5 impossible, use 10/10 = 100% -> core)
        // Actually use another core gene
        for g in 0..10 { matrix.set(g, 1, true); }
        // Shell (3/10 = 30%)
        for g in 0..3 { matrix.set(g, 2, true); }
        // Cloud (1/10 = 10%)
        matrix.set(0, 3, true);

        let temp = NamedTempFile::new().unwrap();
        write_summary_stats(&matrix, temp.path()).unwrap();

        let content = std::fs::read_to_string(temp.path()).unwrap();
        assert!(content.contains("Core genes (>= 99%): 2"));
        assert!(content.contains("Shell genes (15-95%): 1"));
        assert!(content.contains("Cloud genes (0-15%): 1"));
        assert!(content.contains("Total genes: 4"));
    }

    #[test]
    fn test_write_summary_stats_empty() {
        let matrix = BitPackedMatrix::new(0, 0);
        let temp = NamedTempFile::new().unwrap();
        write_summary_stats(&matrix, temp.path()).unwrap();
        // Should succeed without writing (0 genomes)
    }
}