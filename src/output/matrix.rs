//! Presence/absence matrix output (CSV/TSV).

use std::path::Path;
use csv::Writer;

use crate::error::Result;
use crate::graph::BitPackedMatrix;

/// Writer for presence/absence matrix output.
pub struct MatrixWriter;

impl MatrixWriter {
    /// Write the presence/absence matrix to CSV.
    pub fn write(matrix: &BitPackedMatrix, path: &Path) -> Result<()> {
        matrix.to_csv(path)
    }

    /// Write the presence/absence matrix to TSV.
    pub fn write_tsv(matrix: &BitPackedMatrix, path: &Path) -> Result<()> {
        matrix.to_tsv(path)
    }

    /// Write Roary/Panaroo-compatible gene_presence_absence.csv
    ///
    /// Format: Gene, Non-unique Gene name, Annotation, No. isolates,
    /// No. sequences, Avg sequences per isolate, Genome Fragment,
    /// Order within Fragment, Accessory Fragment, Accessory Order with Fragment,
    /// QC, Min group size nuc, Max group size nuc, Avg group size nuc,
    /// then per-isolate columns with gene names or empty.
    pub fn write_roary_csv(matrix: &BitPackedMatrix, path: &Path) -> Result<()> {
        let mut writer = Writer::from_path(path)?;

        // Header: 14 metadata columns + per-genome columns
        let mut header = vec![
            "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
            "No. sequences", "Avg sequences per isolate", "Genome Fragment",
            "Order within Fragment", "Accessory Fragment",
            "Accessory Order with Fragment", "QC",
            "Min group size nuc", "Max group size nuc", "Avg group size nuc",
        ];
        header.extend(matrix.genome_names.iter().map(|s| s.as_str()));

        writer.write_record(&header)?;

        // Data rows
        for (cluster_idx, cluster_id) in matrix.cluster_ids.iter().enumerate() {
            let num_present = matrix.count_present(cluster_idx);
            let avg = if num_present > 0 {
                format!("{:.2}", 1.0) // each genome has exactly 1 copy
            } else {
                "0.00".to_string()
            };

            let mut row = vec![
                cluster_id.clone(),           // Gene
                cluster_id.clone(),           // Non-unique Gene name (same as Gene for now)
                String::new(),                // Annotation (empty)
                num_present.to_string(),      // No. isolates
                num_present.to_string(),      // No. sequences
                avg,                          // Avg sequences per isolate
                String::new(),                // Genome Fragment
                String::new(),                // Order within Fragment
                String::new(),                // Accessory Fragment
                String::new(),                // Accessory Order with Fragment
                String::new(),                // QC
                String::new(),                // Min group size nuc
                String::new(),                // Max group size nuc
                String::new(),                // Avg group size nuc
            ];

            for genome_idx in 0..matrix.num_genomes() {
                row.push(if matrix.get(genome_idx, cluster_idx) {
                    cluster_id.clone()
                } else {
                    String::new()
                });
            }

            writer.write_record(&row)?;
        }

        writer.flush()?;
        Ok(())
    }
}
