//! Panaroo-style output for compatibility and programmatic access.

use std::io::Write;
use std::path::Path;
use serde::Serialize;

use crate::error::Result;
use crate::graph::{PangenomeGraph, BitPackedMatrix, StructuralVariantDetector};

/// Writer for JSON output.
pub struct JsonWriter;

/// JSON representation of the pangenome summary.
#[derive(Serialize)]
struct PangenomeSummary {
    version: String,
    num_genomes: usize,
    num_clusters: usize,
    num_core: usize,
    num_accessory: usize,
    num_edges: usize,
    clusters: Vec<ClusterSummary>,
}

/// JSON representation of a cluster.
#[derive(Serialize)]
struct ClusterSummary {
    id: String,
    support: usize,
    is_paralog: bool,
    is_highly_variable: bool,
    annotations: Vec<String>,
    neighbors: Vec<String>,
}

impl JsonWriter {
    /// Write the pangenome graph and matrix to JSON.
    pub fn write(
        graph: &PangenomeGraph,
        _matrix: &BitPackedMatrix,
        path: &Path,
    ) -> Result<()> {
        let total_genomes = graph.genomes.len().max(1);
        let core_threshold = (total_genomes as f64 * 0.99).ceil() as usize;

        let mut clusters: Vec<ClusterSummary> = Vec::new();
        let mut num_core = 0;

        for (cluster_id, node) in &graph.nodes {
            if node.support >= core_threshold {
                num_core += 1;
            }

            let neighbors: Vec<String> = graph.neighbors(cluster_id)
                .iter()
                .map(|id| id.to_string())
                .collect();

            clusters.push(ClusterSummary {
                id: cluster_id.to_string(),
                support: node.support,
                is_paralog: node.is_paralog,
                is_highly_variable: node.is_highly_variable,
                annotations: node.annotations.iter().cloned().collect(),
                neighbors,
            });
        }

        // Sort clusters by support (descending)
        clusters.sort_by(|a, b| b.support.cmp(&a.support));

        let summary = PangenomeSummary {
            version: crate::VERSION.to_string(),
            num_genomes: total_genomes,
            num_clusters: graph.node_count(),
            num_core,
            num_accessory: graph.node_count() - num_core,
            num_edges: graph.edge_count(),
            clusters,
        };

        let json = serde_json::to_string_pretty(&summary)?;
        std::fs::write(path, json)?;

        Ok(())
    }

    /// Write the gene data CSV (Panaroo-style).
    ///
    /// Links gene sequences to their annotations, with DNA/protein sequences
    /// and location info populated from gene_lookup when available.
    pub fn write_gene_data(
        graph: &PangenomeGraph,
        gene_lookup: &std::collections::HashMap<crate::graph::GeneId, crate::graph::Gene>,
        path: &Path,
    ) -> Result<()> {
        use std::io::Write;
        let mut file = std::fs::File::create(path)?;
        let mut writer = std::io::BufWriter::new(&mut file);

        // Header with DNA and protein sequences
        writeln!(writer, "gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence")?;

        for (_id, node) in &graph.nodes {
            let annotation = node.annotations.iter().next()
                .map(|s| s.as_str())
                .unwrap_or("hypothetical protein");

            // Get location info from first gene member in gene_lookup
            let (contig, start, end, strand): (String, String, String, String) =
                node.gene_members.values().flatten()
                    .filter_map(|gid| gene_lookup.get(&crate::graph::GeneId::new(gid)))
                    .next()
                    .map(|g| (g.contig.clone(), g.start.to_string(), g.end.to_string(), format!("{}", g.strand)))
                    .unwrap_or(("NA".to_string(), "NA".to_string(), "NA".to_string(), "NA".to_string()));

            let dna_seq = node.centroid_sequences.first()
                .map(|s| String::from_utf8_lossy(s).to_string())
                .unwrap_or_default();
            let protein_seq = if let Some(seq) = node.centroid_sequences.first() {
                let protein = crate::io::translate(seq);
                String::from_utf8_lossy(&protein).to_string()
            } else {
                String::new()
            };

            writeln!(writer, "{},{},{},{},{},{},{},{},{},{}",
                node.cluster_id, "", annotation,
                contig, start, end, strand,
                node.support, dna_seq, protein_seq)?;
        }

        Ok(())
    }

    /// Write the pan-genome reference FASTA (Panaroo-style).
    ///
    /// Writes the centroid DNA sequence per cluster.
    /// Clusters without centroid sequences are skipped.
    pub fn write_reference(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        for (cluster_id, node) in &graph.nodes {
            let annotation = node.annotations.iter().next().cloned()
                .unwrap_or_else(|| "hypothetical protein".to_string());

            if let Some(seq) = node.centroid_sequences.first() {
                writeln!(file, ">{} {}", cluster_id, annotation)?;
                writeln!(file, "{}", String::from_utf8_lossy(seq))?;
            }
        }

        Ok(())
    }

    /// Write combined DNA CDS FASTA (Panaroo-style).
    ///
    /// Writes the centroid DNA sequence for each cluster.
    pub fn write_dna_fasta(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        for (cluster_id, node) in &graph.nodes {
            let annotation = node.annotations.iter().next().cloned()
                .unwrap_or_else(|| "hypothetical protein".to_string());

            if let Some(seq) = node.centroid_sequences.first() {
                writeln!(file, ">{} {}", cluster_id, annotation)?;
                writeln!(file, "{}", String::from_utf8_lossy(seq))?;
            }
        }

        Ok(())
    }

    /// Write combined protein CDS FASTA (Panaroo-style).
    ///
    /// Translates the centroid DNA sequence to protein for each cluster.
    pub fn write_protein_fasta(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        for (cluster_id, node) in &graph.nodes {
            let annotation = node.annotations.iter().next().cloned()
                .unwrap_or_else(|| "hypothetical protein".to_string());

            if let Some(seq) = node.centroid_sequences.first() {
                let protein = crate::io::translate(seq);
                if !protein.is_empty() {
                    writeln!(file, ">{} {}", cluster_id, annotation)?;
                    writeln!(file, "{}", String::from_utf8_lossy(&protein))?;
                }
            }
        }

        Ok(())
    }

    /// Write JSON output in compact format (for programmatic access).
    pub fn write_json(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        use serde::Serialize;

        #[derive(Serialize)]
        struct Summary {
            version: String,
            num_genomes: usize,
            num_clusters: usize,
            num_core: usize,
        }

        let total_genomes = graph.genomes.len().max(1);
        let core_threshold = (total_genomes as f64 * 0.99).ceil() as usize;

        let mut num_core = 0;
        for node in graph.nodes.values() {
            if node.support >= core_threshold {
                num_core += 1;
            }
        }

        let summary = Summary {
            version: crate::VERSION.to_string(),
            num_genomes: total_genomes,
            num_clusters: graph.node_count(),
            num_core,
        };

        let json = serde_json::to_string_pretty(&summary)?;
        std::fs::write(path, json)?;

        Ok(())
    }

    /// Write structural variant matrix.
    pub fn write_structural_variants(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        // Use the structural variant detector to find variants
        let detector = StructuralVariantDetector::new();
        let variants = detector.detect(graph);

        // Write the variants to CSV
        crate::output::struct_csv::write_structural_variants(&variants, path)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Node, GeneCluster, GeneId, Gene, GenomeId};
    use tempfile::NamedTempFile;
    use std::collections::HashMap;

    #[test]
    fn test_json_output() {
        let mut graph = PangenomeGraph::new();

        let node = Node::from_cluster(&{
            let mut c = GeneCluster::new("c1");
            c.support = 5;
            c
        });
        graph.add_node(node);

        let matrix = BitPackedMatrix::new(5, 1);

        let temp = NamedTempFile::new().unwrap();
        JsonWriter::write(&graph, &matrix, temp.path()).unwrap();

        let content = std::fs::read_to_string(temp.path()).unwrap();
        assert!(content.contains("num_clusters"));
        assert!(content.contains("c1"));
    }

    #[test]
    fn test_write_gene_data_with_sequences() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("gene_data.csv");

        let mut graph = PangenomeGraph::new();
        let mut node = Node::from_cluster(&GeneCluster::new("c1"));
        node.centroid_sequences = vec![b"ATGCGT".to_vec()];
        node.annotations.insert("hypothetical protein".to_string());
        let mut gene = Gene::new("geneA", GenomeId::new("genome1"));
        gene.contig = "contig1".to_string();
        gene.start = 100;
        gene.end = 105;
        node.gene_members.insert(GenomeId::new("genome1"), vec!["geneA".to_string()]);
        graph.add_node(node);

        let mut gene_lookup = HashMap::new();
        gene_lookup.insert(GeneId::new("geneA"), gene);

        JsonWriter::write_gene_data(&graph, &gene_lookup, &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("dna_sequence"), "Header should have dna_sequence column");
        assert!(content.contains("protein_sequence"), "Header should have protein_sequence column");
        assert!(content.contains("ATGCGT"), "Should contain the DNA sequence");
        assert!(content.contains("contig1"), "Should contain actual contig name");
        assert!(content.contains("100"), "Should contain actual start position");
    }
}
