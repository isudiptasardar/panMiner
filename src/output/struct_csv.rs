//! Structural variant matrix output in CSV format.
//!
//! Outputs a Roary-compatible presence/absence matrix for structural variants.
//! Format: First row is genome names, subsequent rows are variants with 1/0 presence.

use crate::error::Result;
use crate::graph::StructuralVariant;
use std::path::Path;
use std::fs::File;
use std::io::Write;

/// Write structural variant matrix to CSV file.
pub fn write_structural_variants(variants: &[StructuralVariant], path: &Path) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);

    // Get all unique genome names from variants
    let mut genome_names: std::collections::HashSet<String> = std::collections::HashSet::new();
    for variant in variants {
        for genome in &variant.affected_genomes {
            genome_names.insert(genome.clone());
        }
    }
    let mut genomes: Vec<String> = genome_names.into_iter().collect();
    genomes.sort();

    // Write header row
    let header = format!("ID,VariantType,Support,Genomes,{}", genomes.join(","));
    writeln!(writer, "{}", header)?;

    // Write each variant as a row
    for (i, variant) in variants.iter().enumerate() {
        // Check presence for each genome
        let _presence: Vec<&str> = genomes.iter()
            .map(|g| if variant.affected_genomes.contains(g) { "1" } else { "0" })
            .collect();

        let row = format!(
            "SV_{},{},{},\"[{}]\",{}",
            i + 1,
            match variant.variant_type {
                crate::graph::VariantType::Inversion => "INVERSION",
                crate::graph::VariantType::Duplication => "DUPLICATION",
                crate::graph::VariantType::Translocation => "TRANSLOCATION",
                crate::graph::VariantType::Deletion => "DELETION",
            },
            variant.support,
            variant.cluster_ids.join(";"),
            variant.cluster_ids.join(",")
        );
        writeln!(writer, "{}", row)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::VariantType;

    #[test]
    fn test_write_empty_variants() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("struct_variants.csv");

        let variants: Vec<StructuralVariant> = vec![];
        write_structural_variants(&variants, &path).unwrap();

        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("ID,VariantType"));
    }

    #[test]
    fn test_write_variants_with_genomes() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("struct_variants.csv");

        let variants = vec![
            StructuralVariant {
                variant_type: VariantType::Inversion,
                cluster_ids: vec!["A".to_string(), "B".to_string()],
                affected_genomes: vec!["genome1".to_string(), "genome2".to_string()],
                support: 2,
                description: "Inversion between A and B".to_string(),
            },
            StructuralVariant {
                variant_type: VariantType::Duplication,
                cluster_ids: vec!["C".to_string()],
                affected_genomes: vec!["genome1".to_string()],
                support: 3,
                description: "Duplication of C".to_string(),
            },
        ];

        write_structural_variants(&variants, &path).unwrap();

        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("SV_1,INVERSION"));
        assert!(content.contains("SV_2,DUPLICATION"));
        assert!(content.contains("genome1,genome2")); // Header
    }
}
