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
        let core_threshold = (total_genomes as f32 * 0.99).ceil() as usize;

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
    /// Links gene sequences to their annotations.
    pub fn write_gene_data(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        // Write header
        writeln!(file, "gene_id,gene_name,annotation,contig,start,end,strand,support")?;

        // Write each gene in each cluster
        for (cluster_id, node) in &graph.nodes {
            let annotation = node.annotations.iter().next().cloned().unwrap_or_else(|| "hypothetical protein".to_string());
            // Write node info (representative gene for the cluster)
            writeln!(file, "{},,{},NA,NA,NA,NA,{}", cluster_id, annotation, node.support)?;
        }

        Ok(())
    }

    /// Write the pan-genome reference FASTA (Panaroo-style).
    ///
    /// Linear reference genome of all genes found. Paralogous clusters
    /// represented only once to avoid multi-mapping issues.
    pub fn write_reference(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        // Note: This is a placeholder that writes metadata.
        // In a full implementation, this would write actual sequences from the clusters.
        for (cluster_id, node) in &graph.nodes {
            let annotation = node.annotations.iter().next().cloned().unwrap_or_else(|| "hypothetical protein".to_string());
            writeln!(file, ">{} {}", cluster_id, annotation)?;
            writeln!(file, "N{}", "A".repeat(100))?; // Placeholder sequence
        }

        Ok(())
    }

    /// Write combined DNA CDS FASTA (Panaroo-style).
    pub fn write_dna_fasta(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        for (cluster_id, node) in &graph.nodes {
            let annotation = node.annotations.iter().next().cloned().unwrap_or_else(|| "hypothetical protein".to_string());
            writeln!(file, ">{} {}", cluster_id, annotation)?;
            writeln!(file, "N{}", "A".repeat(100))?; // Placeholder sequence
        }

        Ok(())
    }

    /// Write combined protein CDS FASTA (Panaroo-style).
    pub fn write_protein_fasta(graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;

        for (cluster_id, node) in &graph.nodes {
            let annotation = node.annotations.iter().next().cloned().unwrap_or_else(|| "hypothetical protein".to_string());
            writeln!(file, ">{} {}", cluster_id, annotation)?;
            writeln!(file, "M{}", "A".repeat(100))?; // Placeholder sequence starting with Methionine
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
        let core_threshold = (total_genomes as f32 * 0.99).ceil() as usize;

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
    use crate::graph::{Node, GeneCluster};
    use tempfile::NamedTempFile;

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
}
