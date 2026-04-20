//! Extract individual gene sequences from a pangenome output directory.
//!
//! Reads `gene_data.csv` from a PanMiner output directory, filters rows by
//! cluster ID, and writes matching sequences as FASTA.

use std::io::Write;
use std::path::Path;

use crate::error::Result;

/// Row representation of gene_data.csv.
///
/// The current gene_data.csv format has these columns:
/// `gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence`
///
/// Note: in the current output, `gene_id` stores the cluster ID (each row is
/// a cluster centroid), and `gene_name` is typically empty.
#[derive(Debug)]
struct GeneDataRow {
    gene_id: String,
    gene_name: String,
    annotation: String,
    contig: String,
    start: String,
    end: String,
    strand: String,
    support: String,
    dna_sequence: String,
    protein_sequence: String,
}

/// Extract sequences for a given cluster from PanMiner output.
///
/// Reads `gene_data.csv` from the output directory, filters rows where
/// the gene_id column (which stores the cluster ID) matches the given
/// `cluster_id`, and writes matching sequences as FASTA.
///
/// When `protein` is true, writes the protein_sequence column; otherwise
/// writes the dna_sequence column.
pub fn extract_gene(
    output_dir: &Path,
    cluster_id: &str,
    output_path: &Path,
    protein: bool,
) -> Result<()> {
    let gene_data_path = output_dir.join("gene_data.csv");

    if !gene_data_path.exists() {
        return Err(crate::error::Error::InvalidInput(format!(
            "gene_data.csv not found in {:?}",
            output_dir
        )));
    }

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(&gene_data_path)?;

    let headers = reader.headers()?.clone();
    let column_map = build_column_map(&headers);

    let mut file = std::fs::File::create(output_path)?;
    let mut writer = std::io::BufWriter::new(&mut file);
    let mut count = 0usize;

    for result in reader.records() {
        let record = result?;
        let row = parse_row(&record, &column_map)?;

        if row.gene_id != cluster_id {
            continue;
        }

        let sequence = if protein {
            &row.protein_sequence
        } else {
            &row.dna_sequence
        };

        if sequence.is_empty() {
            continue;
        }

        writeln!(
            writer,
            ">{} contig={}",
            row.gene_id, row.contig
        )?;
        writeln!(writer, "{}", sequence)?;
        count += 1;
    }

    if count == 0 {
        return Err(crate::error::Error::ClusterNotFound(cluster_id.to_string()));
    }

    tracing::info!("Extracted {} sequence(s) for cluster '{}'", count, cluster_id);
    Ok(())
}

/// Build a mapping from column name to column index.
fn build_column_map(headers: &csv::StringRecord) -> std::collections::HashMap<String, usize> {
    headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect()
}

/// Get a field from a CSV record by column name, falling back to empty string.
fn get_field<'a>(
    record: &'a csv::StringRecord,
    column_map: &std::collections::HashMap<String, usize>,
    name: &str,
) -> &'a str {
    column_map
        .get(name)
        .and_then(|&i| record.get(i))
        .unwrap_or("")
}

/// Parse a CSV record into a GeneDataRow using the column map.
fn parse_row(
    record: &csv::StringRecord,
    column_map: &std::collections::HashMap<String, usize>,
) -> Result<GeneDataRow> {
    Ok(GeneDataRow {
        gene_id: get_field(record, column_map, "gene_id").to_string(),
        gene_name: get_field(record, column_map, "gene_name").to_string(),
        annotation: get_field(record, column_map, "annotation").to_string(),
        contig: get_field(record, column_map, "contig").to_string(),
        start: get_field(record, column_map, "start").to_string(),
        end: get_field(record, column_map, "end").to_string(),
        strand: get_field(record, column_map, "strand").to_string(),
        support: get_field(record, column_map, "support").to_string(),
        dna_sequence: get_field(record, column_map, "dna_sequence").to_string(),
        protein_sequence: get_field(record, column_map, "protein_sequence").to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn write_gene_data_csv(dir: &Path, content: &str) {
        let path = dir.join("gene_data.csv");
        std::fs::write(path, content).unwrap();
    }

    #[test]
    fn test_extract_gene_finds_cluster() {
        let dir = TempDir::new().unwrap();
        write_gene_data_csv(
            dir.path(),
            "gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence\n\
             cluster_001,,hypothetical protein,contig1,100,200,+,5,ATGCGT,MVR\n\
             cluster_002,,transposase,contig2,300,400,-,2,GGCCAA,GAK\n\
             cluster_003,,hypothetical protein,contig3,500,600,+,1,TTTTTT,FFF\n",
        );

        let output_path = dir.path().join("extracted.fasta");
        extract_gene(dir.path(), "cluster_001", &output_path, false).unwrap();

        let content = std::fs::read_to_string(&output_path).unwrap();
        assert!(content.contains(">cluster_001"));
        assert!(content.contains("ATGCGT"));
        assert!(!content.contains("GGCCAA"));
    }

    #[test]
    fn test_extract_gene_protein_mode() {
        let dir = TempDir::new().unwrap();
        write_gene_data_csv(
            dir.path(),
            "gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence\n\
             cluster_001,,hypothetical protein,contig1,100,200,+,5,ATGCGT,MVR\n",
        );

        let output_path = dir.path().join("extracted.fasta");
        extract_gene(dir.path(), "cluster_001", &output_path, true).unwrap();

        let content = std::fs::read_to_string(&output_path).unwrap();
        assert!(content.contains(">cluster_001"));
        assert!(content.contains("MVR"));
        assert!(!content.contains("ATGCGT"));
    }

    #[test]
    fn test_extract_gene_cluster_not_found() {
        let dir = TempDir::new().unwrap();
        write_gene_data_csv(
            dir.path(),
            "gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence\n\
             cluster_001,,hypothetical protein,contig1,100,200,+,5,ATGCGT,MVR\n",
        );

        let output_path = dir.path().join("extracted.fasta");
        let result = extract_gene(dir.path(), "cluster_999", &output_path, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_extract_gene_missing_file() {
        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("extracted.fasta");
        let result = extract_gene(dir.path(), "cluster_001", &output_path, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_extract_gene_skips_empty_sequences() {
        let dir = TempDir::new().unwrap();
        write_gene_data_csv(
            dir.path(),
            "gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence\n\
             cluster_001,,hypothetical protein,contig1,100,200,+,5,,MVR\n\
             cluster_001,,hypothetical protein,contig2,300,400,+,5,ATGCGT,MVR\n",
        );

        let output_path = dir.path().join("extracted.fasta");
        extract_gene(dir.path(), "cluster_001", &output_path, false).unwrap();

        let content = std::fs::read_to_string(&output_path).unwrap();
        assert!(content.contains("ATGCGT"));
        // The empty-DNA row should not produce a FASTA entry in DNA mode
        assert_eq!(content.matches(">").count(), 1);
    }
}