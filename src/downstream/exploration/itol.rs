//! iTOL annotation export for phylogenetic tree visualization.
//!
//! Exports phylogenetic tree annotations for iTOL (Interactive Tree of Life).
//! Generates color strip, text label, and binary presence heatmap datasets.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};
use crate::error::Result;
use crate::graph::{ClusterId, GenomeId, Node, PangenomeGraph};

/// Export phylogenetic tree annotations for iTOL visualization.
#[derive(Debug, Clone)]
pub struct ItolAnnotationRunner {
    tree_file: Option<PathBuf>,
    color_by_source: bool,
}

impl ItolAnnotationRunner {
    /// Create a new iTOL annotation runner.
    pub fn new() -> Self {
        Self {
            tree_file: None,
            color_by_source: false,
        }
    }

    /// Set the phylogenetic tree file (Newick format).
    pub fn with_tree(mut self, path: PathBuf) -> Self {
        self.tree_file = Some(path);
        self
    }

    /// Enable coloring by genome source file.
    pub fn with_color_by_source(mut self, enabled: bool) -> Self {
        self.color_by_source = enabled;
        self
    }

    /// Parse GML file into a PangenomeGraph.
    fn parse_gml_graph(gml_path: &Path) -> Result<PangenomeGraph> {
        let content = fs::read_to_string(gml_path)?;
        let mut graph = PangenomeGraph::new();

        let mut current_node: Option<Node> = None;
        let mut in_node = false;
        let mut in_genomes = false;

        for line in content.lines() {
            let line = line.trim();

            if line == "node [" {
                in_node = true;
                current_node = None;
            } else if line == "]" && in_genomes {
                in_genomes = false;
            } else if line == "]" && in_node {
                if let Some(node) = current_node.take() {
                    graph.nodes.insert(node.cluster_id.clone(), node);
                }
                in_node = false;
            } else if in_node {
                if line.starts_with("label") {
                    let label = line.split('"').nth(1).unwrap_or("").to_string();
                    let mut node = Node::from_cluster(&{
                        let c = crate::graph::GeneCluster::new(&label);
                        c
                    });
                    node.cluster_id = ClusterId::new(&label);
                    current_node = Some(node);
                } else if line.starts_with("support") {
                    if let Some(s) = line.split_whitespace().nth(1) {
                        if let Some(node) = &mut current_node {
                            node.support = s.parse().unwrap_or(1);
                        }
                    }
                } else if line.starts_with("genomes [") {
                    in_genomes = true;
                } else if in_genomes && line.starts_with('"') {
                    let genome_id = line.trim_matches('"').to_string();
                    if let Some(node) = &mut current_node {
                        node.genomes.insert(GenomeId::new(genome_id));
                    }
                } else if line.starts_with("is_paralog") {
                    if let Some(node) = &mut current_node {
                        if let Some(v) = line.split_whitespace().nth(1) {
                            node.is_paralog = v == "1";
                        }
                    }
                } else if line.starts_with("is_highly_variable") {
                    if let Some(node) = &mut current_node {
                        if let Some(v) = line.split_whitespace().nth(1) {
                            node.is_highly_variable = v == "1";
                        }
                    }
                }
            }
        }

        Ok(graph)
    }

    /// Parse gene_presence_absence.csv (Roary format) to extract per-genome metadata.
    ///
    /// Returns a HashMap mapping genome name to gene count string.
    /// Skips the 14 metadata columns in the Roary CSV header.
    #[allow(dead_code)]
    fn parse_gene_data_csv(csv_path: &Path) -> Result<HashMap<String, String>> {
        // Expects path to gene_presence_absence.csv (Roary format: 14 metadata cols + genome cols)
        let content = fs::read_to_string(csv_path)?;
        let mut metadata: HashMap<String, String> = HashMap::new();

        let mut lines = content.lines();
        if let Some(header) = lines.next() {
            let headers: Vec<&str> = header.split(',').collect();
            const METADATA_COLS: usize = 14;
            if headers.len() <= METADATA_COLS {
                return Ok(metadata); // No genome columns
            }

            let genome_names: Vec<&str> = headers[METADATA_COLS..]
                .iter()
                .map(|s| s.trim())
                .collect();

            // Initialize gene counts per genome
            let mut gene_counts: Vec<usize> = vec![0; genome_names.len()];

            for line in lines {
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                let fields: Vec<&str> = line.split(',').collect();
                for (idx, _name) in genome_names.iter().enumerate() {
                    if let Some(val) = fields.get(METADATA_COLS + idx) {
                        let val = val.trim();
                        if !val.is_empty() && val != "0" {
                            gene_counts[idx] += 1;
                        }
                    }
                }
            }

            for (name, count) in genome_names.iter().zip(gene_counts.iter()) {
                metadata.insert(name.to_string(), count.to_string());
            }
        }

        Ok(metadata)
    }

    /// Build an allelic profile matrix from the graph.
    fn build_profile_matrix(
        graph: &PangenomeGraph,
    ) -> (Vec<String>, Vec<String>, Vec<Vec<u8>>) {
        let mut genome_ids: Vec<String> = graph
            .nodes
            .values()
            .flat_map(|node| node.genomes.iter())
            .map(|g| g.as_str().to_string())
            .collect();
        genome_ids.sort();
        genome_ids.dedup();

        let mut cluster_ids: Vec<String> = graph
            .nodes
            .keys()
            .map(|c| c.as_str().to_string())
            .collect();
        cluster_ids.sort();

        let mut presence = vec![vec![0u8; cluster_ids.len()]; genome_ids.len()];

        for (cluster_idx, cluster_id) in cluster_ids.iter().enumerate() {
            let cluster_id = ClusterId::new(cluster_id);
            if let Some(node) = graph.nodes.get(&cluster_id) {
                for (genome_idx, genome_id) in genome_ids.iter().enumerate() {
                    let genome_id = GenomeId::new(genome_id);
                    if node.genomes.contains(&genome_id) {
                        presence[genome_idx][cluster_idx] = 1;
                    }
                }
            }
        }

        (genome_ids, cluster_ids, presence)
    }

    /// Write the tree file (copy Newick to output).
    fn write_itol_tree(tree_path: &Path, output_dir: &Path) -> Result<PathBuf> {
        if !tree_path.exists() {
            return Err(crate::error::Error::InvalidInput(format!(
                "Tree file not found: {}",
                tree_path.display()
            )));
        }
        let dest_path = output_dir.join("itol_tree.txt");
        fs::copy(tree_path, &dest_path)?;
        Ok(dest_path)
    }

    /// Write iTOL annotations with color strips and binary presence heatmap.
    fn write_itol_annotations(
        genome_ids: &[String],
        cluster_ids: &[String],
        presence: &[Vec<u8>],
        output_dir: &Path,
    ) -> Result<PathBuf> {
        let itol_path = output_dir.join("itol_annotations.txt");
        let mut file = fs::File::create(&itol_path)?;
        use std::io::Write;

        let num_genomes = genome_ids.len();

        // Determine core vs accessory genes
        let mut core_clusters = vec![false; cluster_ids.len()];
        for cluster_idx in 0..cluster_ids.len() {
            let all_present = presence.iter().all(|row| row[cluster_idx] == 1);
            core_clusters[cluster_idx] = all_present;
        }

        // --- DATASET STIPED: Core vs Accessory color strip ---
        writeln!(file, "DATASET STIPED")?;
        writeln!(file, "# Color strip dataset for gene categories")?;
        writeln!(file, "DATASET_TITLE\tCore vs Accessory Genes")?;
        writeln!(file, "COLOR\t#ff0000")?;
        writeln!(file)?;

        writeln!(file, "LEGEND_TITLE\tGene Category")?;
        writeln!(file, "LEGEND_COLORS\t#ff0000\t#008080\t#FFA500")?;
        writeln!(file, "LEGEND_LABELS\tCore\tAccessory\tParalog")?;
        writeln!(file, "DATA")?;

        for (genome_idx, genome_id) in genome_ids.iter().enumerate() {
            let mut core_count = 0;
            let mut acc_count = 0;
            for cluster_idx in 0..cluster_ids.len() {
                if presence[genome_idx][cluster_idx] == 1 {
                    if core_clusters[cluster_idx] {
                        core_count += 1;
                    } else {
                        acc_count += 1;
                    }
                }
            }
            let color = if core_count >= acc_count {
                "#ff0000" // core
            } else {
                "#008080" // accessory
            };
            writeln!(file, "{}\t{}", genome_id, color)?;
        }

        writeln!(file)?;

        // --- DATASET SIMPLEBAR: Gene count per genome ---
        writeln!(file, "DATASET SIMPLEBAR")?;
        writeln!(file, "# Gene count per genome")?;
        writeln!(file, "DATASET_TITLE\tGene Count")?;
        writeln!(file, "COLOR\t#008080")?;
        writeln!(file)?;

        writeln!(file, "DATA")?;
        for (genome_idx, genome_id) in genome_ids.iter().enumerate() {
            let gene_count: usize = presence[genome_idx].iter().map(|&v| v as usize).sum();
            writeln!(file, "{}\t{}", genome_id, gene_count)?;
        }

        writeln!(file)?;

        // --- DATASET BINARY: Top N most variable genes heatmap ---
        writeln!(file, "DATASET BINARY")?;
        writeln!(file, "# Binary presence for top variable genes")?;
        writeln!(file, "DATASET_TITLE\tTop Variable Genes")?;
        writeln!(file, "COLOR\t#008080")?;
        writeln!(file)?;

        // Compute variance for each cluster
        let mut cluster_variance: Vec<(usize, f64)> = cluster_ids
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let count: usize = presence.iter().map(|row| row[idx] as usize).sum();
                let p = count as f64 / num_genomes as f64;
                let variance = p * (1.0 - p);
                (idx, variance)
            })
            .collect();

        cluster_variance.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let top_indices: Vec<usize> = cluster_variance
            .into_iter()
            .take(50)
            .map(|(idx, _)| idx)
            .collect();

        // Field labels
        writeln!(file, "DATA")?;
        for (genome_idx, genome_id) in genome_ids.iter().enumerate() {
            let values: Vec<String> = top_indices
                .iter()
                .map(|&idx| presence[genome_idx][idx].to_string())
                .collect();
            writeln!(file, "{}\t{}", genome_id, values.join("|"))?;
        }

        Ok(itol_path)
    }
}

impl Default for ItolAnnotationRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl DownstreamRunner for ItolAnnotationRunner {
    type Output = ItolResult;

    fn run(&self, output_dir: &Path) -> Result<ItolResult> {
        let gml_path = output_dir.join("final_graph.gml");
        if !gml_path.exists() {
            return Err(crate::error::Error::InvalidInput(format!(
                "final_graph.gml not found at {:?}",
                gml_path
            )));
        }

        let graph = Self::parse_gml_graph(&gml_path)?;
        if graph.nodes.is_empty() {
            return Err(crate::error::Error::InvalidInput(
                "No nodes found in graph".to_string(),
            ));
        }

        let (genome_ids, cluster_ids, presence) = Self::build_profile_matrix(&graph);

        let downstream_dir = output_dir.join("downstream");
        fs::create_dir_all(&downstream_dir)?;

        // Write iTOL tree (copy from input if provided)
        let itol_tree_path = if let Some(ref tree_path) = self.tree_file {
            Some(Self::write_itol_tree(tree_path, &downstream_dir)?)
        } else {
            None
        };

        // Write iTOL annotations
        let itol_path = Self::write_itol_annotations(&genome_ids, &cluster_ids, &presence, &downstream_dir)?;

        Ok(ItolResult {
            itol_tree_path,
            itol_annotations_path: itol_path,
            num_genomes: genome_ids.len(),
            num_genes: cluster_ids.len(),
        })
    }

    fn name(&self) -> &str {
        "iTOL Export"
    }

    fn is_available(&self) -> bool {
        true // Always available - pure export
    }

    fn required_inputs(&self) -> Vec<DownstreamInput> {
        vec![DownstreamInput::FinalGraph]
    }
}

/// Result of iTOL annotation export.
#[derive(Debug)]
pub struct ItolResult {
    pub itol_tree_path: Option<PathBuf>,
    pub itol_annotations_path: PathBuf,
    pub num_genomes: usize,
    pub num_genes: usize,
}

impl DownstreamResult for ItolResult {
    fn write_to(&self, _dir: &Path) -> Result<()> {
        if !self.itol_annotations_path.exists() {
            return Err(crate::error::Error::InvalidInput(format!(
                "iTOL annotations file not found: {:?}",
                self.itol_annotations_path
            )));
        }
        if let Some(ref tree_path) = self.itol_tree_path {
            if !tree_path.exists() {
                return Err(crate::error::Error::InvalidInput(format!(
                    "iTOL tree file not found: {:?}",
                    tree_path
                )));
            }
        }
        Ok(())
    }

    fn summary(&self) -> String {
        format!(
            "iTOL export: {} genomes, {} gene clusters, annotations: {}",
            self.num_genomes,
            self.num_genes,
            self.itol_annotations_path.display()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_itol_annotation_result() {
        let temp_dir = TempDir::new().unwrap();
        let result = ItolResult {
            itol_tree_path: None,
            itol_annotations_path: temp_dir.path().join("itol_annotations.txt"),
            num_genomes: 10,
            num_genes: 50,
        };
        // Just verify summary doesn't panic
        let summary = result.summary();
        assert!(summary.contains("iTOL export"));
        assert!(summary.contains("10 genomes"));
    }
}
