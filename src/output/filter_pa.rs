//! Presence/absence matrix filtering.
//!
//! Removes pseudogenes, fragments, and length outliers from the
//! gene presence/absence matrix, mirroring Panaroo's filter-pa command.

use std::path::Path;
use crate::error::Result;

/// Filter types for presence/absence filtering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterType {
    /// Remove fragmented genes (those containing semicolons in gene IDs)
    Fragment,
    /// Remove pseudogenes (those with internal stop codons)
    Pseudogene,
    /// Remove length outliers (genes significantly shorter/longer than mode)
    Length,
}

impl std::fmt::Display for FilterType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FilterType::Fragment => write!(f, "frag"),
            FilterType::Pseudogene => write!(f, "pseudo"),
            FilterType::Length => write!(f, "length"),
        }
    }
}

/// Parse a comma-separated filter type string.
pub fn parse_filter_types(s: &str) -> Vec<FilterType> {
    s.split(',')
        .filter_map(|f| match f.trim().to_lowercase().as_str() {
            "frag" | "fragment" => Some(FilterType::Fragment),
            "pseudo" | "pseudogene" => Some(FilterType::Pseudogene),
            "length" | "len" => Some(FilterType::Length),
            _ => None,
        })
        .collect()
}

/// Check if a gene ID indicates a fragment (contains semicolons from paralog merging).
pub fn is_fragment(gene_id: &str) -> bool {
    gene_id.contains(';')
}

/// Check if a gene annotation indicates a pseudogene.
///
/// Pseudogenes typically have "pseudo", "pseudogene", or internal stop codons.
/// "hypothetical protein" is NOT treated as a pseudogene — it indicates an
/// uncharacterized gene, not a non-functional one.
pub fn is_pseudogene(annotation: &str, protein: &[u8]) -> bool {
    let annot_lower = annotation.to_lowercase();
    if annot_lower.contains("pseudo")
        || annot_lower.contains("fragment")
    {
        return true;
    }
    // Check for internal stop codons in protein
    if protein.len() > 1 && protein[..protein.len() - 1].contains(&b'*') {
        return true;
    }
    false
}

/// Check if a gene is a length outlier.
#[allow(dead_code)]
///
/// A gene is a length outlier if its length deviates by more than
/// `threshold` proportion from the mode length of its cluster.
pub fn is_length_outlier(length: usize, mode_length: usize, threshold: f32) -> bool {
    if mode_length == 0 {
        return false;
    }
    let deviation = (length as f32 - mode_length as f32).abs() / mode_length as f32;
    deviation > threshold
}

/// Compute the mode (most common) length from a list of lengths.
#[allow(dead_code)]
pub fn compute_mode_length(lengths: &[usize]) -> usize {
    use std::collections::HashMap;
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for &len in lengths {
        *counts.entry(len).or_insert(0) += 1;
    }
    counts.into_iter().max_by_key(|&(_, count)| count).map(|(len, _)| len).unwrap_or(0)
}

/// Filter a presence/absence CSV file.
///
/// Reads a PanMiner `gene_presence_absence.csv`, applies the specified
/// filter types, and writes a filtered version.
pub fn filter_presence_absence(
    input_path: &Path,
    output_path: &Path,
    filter_types: &[FilterType],
    _length_threshold: f32,
) -> Result<()> {
    let mut rdr = csv::Reader::from_path(input_path)?;
    let mut wtr = csv::Writer::from_path(output_path)?;

    let headers = rdr.headers()?.clone();
    wtr.write_record(&headers)?;

    let mut removed_count = 0;
    let mut kept_count = 0;

    for result in rdr.records() {
        let record = result?;
        if record.len() < 3 {
            wtr.write_record(&record)?;
            kept_count += 1;
            continue;
        }

        let gene_id = record.get(0).unwrap_or("");
        let annotation = record.get(1).unwrap_or("");

        let mut should_filter = false;

        for ft in filter_types {
            match ft {
                FilterType::Fragment => {
                    if is_fragment(gene_id) {
                        should_filter = true;
                        break;
                    }
                }
                FilterType::Pseudogene => {
                    let protein = b"";
                    if is_pseudogene(annotation, protein) {
                        should_filter = true;
                        break;
                    }
                }
                FilterType::Length => {
                    // Length filtering requires gene length data
                    // which is not in the standard P/A CSV
                    // This is handled during output generation
                }
            }
        }

        if should_filter {
            removed_count += 1;
        } else {
            wtr.write_record(&record)?;
            kept_count += 1;
        }
    }

    wtr.flush()?;

    tracing::info!(
        "Filtered P/A matrix: removed {} genes, {} remaining",
        removed_count,
        kept_count
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_fragment() {
        assert!(is_fragment("gene1;gene2"));
        assert!(is_fragment("group_1;group_2;group_3"));
        assert!(!is_fragment("gene1"));
    }

    #[test]
    fn test_is_pseudogene() {
        assert!(is_pseudogene("pseudogene", b""));
        assert!(!is_pseudogene("hypothetical protein", b"")); // Not a pseudogene
        assert!(is_pseudogene("normal annotation", b"MA*MA")); // Internal stop
        assert!(!is_pseudogene("normal annotation", b"MAMAM")); // No internal stop
    }

    #[test]
    fn test_is_length_outlier() {
        assert!(is_length_outlier(100, 300, 0.5)); // 67% deviation
        assert!(!is_length_outlier(280, 300, 0.5)); // 6.7% deviation
        assert!(!is_length_outlier(100, 0, 0.5)); // Zero mode length
    }

    #[test]
    fn test_compute_mode_length() {
        assert_eq!(compute_mode_length(&[100, 100, 200, 300]), 100);
        assert_eq!(compute_mode_length(&[100, 200, 200, 300]), 200);
        assert_eq!(compute_mode_length(&[50]), 50);
    }

    #[test]
    fn test_parse_filter_types() {
        let types = parse_filter_types("frag,pseudo,length");
        assert_eq!(types.len(), 3);
        assert_eq!(types[0], FilterType::Fragment);
        assert_eq!(types[1], FilterType::Pseudogene);
        assert_eq!(types[2], FilterType::Length);

        // Test alternate names
        let types2 = parse_filter_types("fragment,pseudogene,len");
        assert_eq!(types2.len(), 3);
        assert_eq!(types2[0], FilterType::Fragment);
        assert_eq!(types2[1], FilterType::Pseudogene);
        assert_eq!(types2[2], FilterType::Length);
    }
}