//! Presence/absence matrix output (CSV/TSV).

use std::collections::HashMap;
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

    /// Write a Roary-compatible gene P/A CSV with semicolon-delimited gene member IDs.
    ///
    /// Same 14-column header as `write_roary_csv`, but per-genome cells contain
    /// semicolon-delimited gene IDs (e.g., "geneA;geneB") instead of just the cluster ID.
    /// This matches Panaroo's `gene_presence_absence_roary.csv` format.
    pub fn write_roary_gene_csv(
        matrix: &BitPackedMatrix,
        gene_members: &HashMap<String, HashMap<String, Vec<String>>>,
        path: &Path,
    ) -> Result<()> {
        let mut wtr = csv::Writer::from_path(path)?;

        // Header: 14 fixed columns + per-genome columns
        let mut header: Vec<String> = vec![
            "Gene".into(), "Non-unique Gene name".into(), "Annotation".into(),
            "No. isolates".into(), "No. sequences".into(),
            "Avg sequences per isolate".into(), "Genome Fragment".into(),
            "Order within Fragment".into(), "Accessory Fragment".into(),
            "Accessory Order with Fragment".into(), "QC".into(),
            "Min group size nuc".into(), "Max group size nuc".into(),
            "Avg group size nuc".into(),
        ];
        header.extend(matrix.genome_names.iter().cloned());
        wtr.write_record(&header)?;

        // Data rows
        for (cluster_idx, cluster_id) in matrix.cluster_ids.iter().enumerate() {
            let count_present = matrix.count_present(cluster_idx);
            let avg = if count_present > 0 { "1.00" } else { "0.00" };

            let mut row: Vec<String> = vec![
                cluster_id.clone(), cluster_id.clone(), String::new(),
                count_present.to_string(), count_present.to_string(), avg.into(),
                String::new(), String::new(), String::new(), String::new(),
                String::new(), String::new(), String::new(), String::new(),
            ];

            let members = gene_members.get(cluster_id);

            for genome_idx in 0..matrix.num_genomes() {
                let genome_name = &matrix.genome_names[genome_idx];
                if matrix.get(genome_idx, cluster_idx) {
                    if let Some(members) = members {
                        if let Some(gene_ids) = members.get(genome_name) {
                            row.push(gene_ids.join(";"));
                        } else {
                            row.push(cluster_id.clone());
                        }
                    } else {
                        row.push(cluster_id.clone());
                    }
                } else {
                    row.push(String::new());
                }
            }
            wtr.write_record(&row)?;
        }

        wtr.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_roary_gene_csv() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("roary_gene.csv");

        let mut matrix = BitPackedMatrix::new(2, 1);
        matrix.set_genome_names(vec!["genome1".to_string(), "genome2".to_string()]);
        matrix.set_cluster_ids(vec!["cluster_0".to_string()]);
        matrix.set(0, 0, true);
        matrix.set(1, 0, true);

        let mut gene_members: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        let mut members = HashMap::new();
        members.insert("genome1".to_string(), vec!["geneA".to_string()]);
        members.insert("genome2".to_string(), vec!["geneB".to_string(), "geneC".to_string()]);
        gene_members.insert("cluster_0".to_string(), members);

        MatrixWriter::write_roary_gene_csv(&matrix, &gene_members, &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("geneA"), "Should contain geneA in genome1 column");
        assert!(content.contains("geneB;geneC"), "Should contain semicolon-delimited genes in genome2 column");
    }
}
