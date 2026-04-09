//! Quality control statistics output.

use crate::error::Result;
use crate::io::GenomeQC;
use std::path::Path;
use std::fs::File;
use std::io::Write;

/// Write QC statistics to CSV file.
pub fn write_qc_stats(qc_results: &[GenomeQC], path: &Path) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);

    // Write header
    writeln!(writer, "GenomeID,Completeness,Contamination,GenomeSize,NumContigs,N50,MashDistance,Passed")?;

    // Write each genome's QC data
    for qc in qc_results {
        let passed = if qc.passed { "true" } else { "false" };
        let mash_dist = qc.mash_distance.map(|d| format!("{:.4}", d)).unwrap_or_else(|| "N/A".to_string());

        writeln!(writer, "{},{:.2},{:.2},{},{},{},{},{}",
            qc.genome_id,
            qc.completeness,
            qc.contamination,
            qc.genome_size,
            qc.num_contigs,
            qc.n50,
            mash_dist,
            passed
        )?;
    }

    Ok(())
}

/// Write QC summary to text file.
pub fn write_qc_summary(qc_results: &[GenomeQC], path: &Path) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);

    let total = qc_results.len();
    let passed: usize = qc_results.iter().filter(|q| q.passed).count();
    let failed = total - passed;

    let avg_completeness = if total > 0 {
        qc_results.iter().map(|q| q.completeness).sum::<f64>() / total as f64
    } else {
        0.0
    };

    let avg_contamination = if total > 0 {
        qc_results.iter().map(|q| q.contamination).sum::<f64>() / total as f64
    } else {
        0.0
    };

    writeln!(writer, "=== PanMiner QC Summary ===")?;
    writeln!(writer, "Total genomes: {}", total)?;
    writeln!(writer, "Passed QC: {}", passed)?;
    writeln!(writer, "Failed QC: {}", failed)?;
    writeln!(writer, "Avg completeness: {:.2}%", avg_completeness)?;
    writeln!(writer, "Avg contamination: {:.2}%", avg_contamination)?;
    writeln!(writer)?;

    writeln!(writer, "=== Genomes ===")?;
    for qc in qc_results {
        writeln!(writer, "{}: {}% complete, {}% contamination - {}",
            qc.genome_id,
            qc.completeness,
            qc.contamination,
            if qc.passed { "PASS" } else { "FAIL" }
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_empty_qc() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("qc.csv");

        let qc_results: Vec<GenomeQC> = vec![];
        write_qc_stats(&qc_results, &path).unwrap();

        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("GenomeID,Completeness,Contamination"));
    }

    #[test]
    fn test_write_qc_summary() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("qc_summary.txt");

        let qc_results = vec![
            GenomeQC {
                genome_id: "genome1".to_string(),
                completeness: 95.0,
                contamination: 2.0,
                passed: true,
                ..Default::default()
            },
            GenomeQC {
                genome_id: "genome2".to_string(),
                completeness: 85.0,
                contamination: 8.0,
                passed: false,
                ..Default::default()
            },
        ];

        write_qc_summary(&qc_results, &path).unwrap();

        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("PanMiner QC Summary"));
        assert!(content.contains("Total genomes: 2"));
    }
}
