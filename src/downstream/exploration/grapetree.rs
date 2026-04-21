//! GrapeTree and iTOL export for pangenome visualization.
//!
//! Exports pangenome data for visualization in GrapeTree (minimum spanning trees)
//! and iTOL (phylogenetic trees with annotation tracks). GrapeTree invocation is
//! optional - profiles are always exported regardless of whether the tool is installed.
//!
//! # Outputs
//!
//! - `grapetree_profiles.tsv` - Allelic profile matrix (GrapeTree input format)
//! - `itol_annotations.txt` - iTOL dataset file (color strips + gene presence heatmap)

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};
use crate::error::{Error, Result};
use crate::graph::{ClusterId, GenomeId, Node, PangenomeGraph};

/// Number of most variable genes to include in iTOL heatmap.
const ITOL_TOP_N_VARIABLE: usize = 50;

/// Export runner for GrapeTree and iTOL visualization.
#[derive(Debug, Clone)]
pub struct GrapeTreeExportRunner {
    /// Include iTOL annotation export.
    include_itol: bool,
    /// Number of top variable genes for iTOL heatmap.
    #[allow(dead_code)]
    top_n_variable: usize,
}

impl GrapeTreeExportRunner {
    /// Create a new GrapeTree export runner.
    pub fn new(include_itol: bool) -> Self {
        Self {
            include_itol,
            top_n_variable: ITOL_TOP_N_VARIABLE,
        }
    }

    /// Detect if GrapeTree is installed and available.
    pub fn detect_grapetree() -> bool {
        which::which("grapetree").is_ok()
    }

    /// Parse the pangenome graph from a GML file.
    fn parse_gml_graph(gml_path: &Path) -> Result<PangenomeGraph> {
        let content = std::fs::read_to_string(gml_path)?;
        let mut graph = PangenomeGraph::new();

        // GML format: node [ id "cluster_id" label "cluster_id" support N genomes [ ... ] ... ]
        //             edge [ source "id" target "id" ... ]
        let mut current_node: Option<Node> = None;
        let mut in_node = false;
        let mut in_edge = false;
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
            } else if line == "edge [" {
                in_edge = true;
            } else if line == "]" && in_edge {
                in_edge = false;
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
                    // Genome ID in quotes
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

    /// Build an allelic profile matrix from the graph.
    ///
    /// Returns: (genome_ids, cluster_ids, presence_matrix)
    ///   where presence_matrix[genome_idx][cluster_idx] is 1 if present, 0 if absent.
    fn build_profile_matrix(
        graph: &PangenomeGraph,
    ) -> (Vec<String>, Vec<String>, Vec<Vec<u8>>) {
        // Collect all genome IDs (sorted for deterministic output)
        let mut genome_ids: Vec<String> = graph
            .nodes
            .values()
            .flat_map(|node| node.genomes.iter())
            .map(|g| g.as_str().to_string())
            .collect();
        genome_ids.sort();
        genome_ids.dedup();

        // Collect all cluster IDs (sorted for deterministic output)
        let mut cluster_ids: Vec<String> = graph
            .nodes
            .keys()
            .map(|c| c.as_str().to_string())
            .collect();
        cluster_ids.sort();

        // Build presence matrix
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

    /// Write GrapeTree profile TSV file.
    fn write_profiles(
        profiles_path: &Path,
        genome_ids: &[String],
        cluster_ids: &[String],
        presence: &[Vec<u8>],
    ) -> Result<()> {
        let mut file = std::fs::File::create(profiles_path)?;
        use std::io::Write;

        // Header: "genome" <tab> cluster1 <tab> cluster2 ...
        write!(file, "genome")?;
        for cluster_id in cluster_ids {
            write!(file, "\t{}", cluster_id)?;
        }
        writeln!(file)?;

        // Rows: genome_id <tab> 0/1 <tab> 0/1 ...
        for (genome_idx, genome_id) in genome_ids.iter().enumerate() {
            write!(file, "{}", genome_id)?;
            for &val in &presence[genome_idx] {
                write!(file, "\t{}", val)?;
            }
            writeln!(file)?;
        }

        Ok(())
    }

    /// Write iTOL annotation file.
    fn write_itol_annotations(
        itol_path: &Path,
        genome_ids: &[String],
        cluster_ids: &[String],
        presence: &[Vec<u8>],
        _total_genomes: usize,
    ) -> Result<()> {
        let mut file = std::fs::File::create(itol_path)?;
        use std::io::Write;

        // --- Color strip dataset: core (red) vs accessory (teal) per genome ---
        writeln!(file, "DATASET STIPED")?;
        writeln!(file, "# Species, strains or any other categories")?;
        writeln!(file, "DATASET_TITLE\tCore vs Accessory Genes")?;
        writeln!(file, "COLOR\t#ff0000")?;
        writeln!(file)?;

        // Determine core genes (present in all genomes)
        let num_genomes = genome_ids.len();
        let mut core_clusters = vec![false; cluster_ids.len()];
        for cluster_idx in 0..cluster_ids.len() {
            let all_present = presence.iter().all(|row| row[cluster_idx] == 1);
            core_clusters[cluster_idx] = all_present;
        }

        // Write legend and data
        writeln!(file, "LEGEND_TITLE\tGene Category")?;
        writeln!(file, "LEGEND_COLORS\t#ff0000\t#008080")?;
        writeln!(file, "LEGEND_LABELS\tCore\tAccessory")?;
        writeln!(file, "DATA")?;

        for (genome_idx, genome_id) in genome_ids.iter().enumerate() {
            // Count core and accessory genes for this genome
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

            // Color: red for core-dominated, teal for accessory-dominated
            let color = if core_count >= acc_count {
                "#ff0000" // core
            } else {
                "#008080" // accessory
            };
            writeln!(file, "{}\t{}", genome_id, color)?;
        }

        writeln!(file)?;

        // --- Simple bar dataset: gene presence for top N variable genes ---
        writeln!(file, "DATASET SIMPLEBAR")?;
        writeln!(file, "# Numeric gene presence values")?;
        writeln!(file, "DATASET_TITLE\tTop Variable Gene Presence")?;
        writeln!(file, "COLOR\t#008080")?;
        writeln!(file)?;

        // Compute variance for each cluster (number of genomes that have it)
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

        // Sort by variance descending and take top N
        cluster_variance.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let top_indices: Vec<usize> = cluster_variance
            .into_iter()
            .take(ITOL_TOP_N_VARIABLE)
            .map(|(idx, _)| idx)
            .collect();

        writeln!(file, "DATA")?;
        for (genome_idx, genome_id) in genome_ids.iter().enumerate() {
            let total_present: f64 = top_indices
                .iter()
                .map(|&idx| presence[genome_idx][idx] as f64)
                .sum();
            writeln!(file, "{}\t{:.1}", genome_id, total_present)?;
        }

        Ok(())
    }

    /// Run GrapeTree if it is installed.
    fn run_grapetree(profiles_path: &Path, output_prefix: &Path) -> Result<()> {
        if !Self::detect_grapetree() {
            return Err(Error::ExternalTool(
                "GrapeTree is not installed. Install with: conda install -c bioconda grapetree".to_string(),
            ));
        }

        let output_prefix_str = output_prefix
            .to_str()
            .ok_or_else(|| Error::Output("Invalid output prefix path".to_string()))?;

        let output = Command::new("grapetree")
            .args(["-i", &profiles_path.to_string_lossy()])
            .args(["-o", output_prefix_str])
            .output()
            .map_err(|e| Error::ExternalTool(format!("Failed to run grapetree: {}", e)))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::ExternalTool(format!(
                "GrapeTree failed: {}",
                stderr
            )));
        }

        Ok(())
    }
}

impl DownstreamRunner for GrapeTreeExportRunner {
    type Output = GrapetreeResult;

    fn run(&self, output_dir: &Path) -> Result<Self::Output> {
        let gml_path = output_dir.join("final_graph.gml");

        if !gml_path.exists() {
            return Err(Error::Output(format!(
                "final_graph.gml not found at {:?}",
                gml_path
            )));
        }

        // Parse the graph
        let graph = Self::parse_gml_graph(&gml_path)?;

        if graph.nodes.is_empty() {
            return Err(Error::Output(
                "No nodes found in graph".to_string(),
            ));
        }

        // Build profile matrix
        let (genome_ids, cluster_ids, presence) = Self::build_profile_matrix(&graph);

        let num_genomes = genome_ids.len();
        let num_genes = cluster_ids.len();

        // Create downstream output directory
        let downstream_dir = output_dir.join("downstream");
        std::fs::create_dir_all(&downstream_dir)?;

        // Write profiles TSV
        let profiles_path = downstream_dir.join("grapetree_profiles.tsv");
        Self::write_profiles(&profiles_path, &genome_ids, &cluster_ids, &presence)?;

        // Optionally write iTOL annotations
        let itol_path = if self.include_itol {
            let path = downstream_dir.join("itol_annotations.txt");
            Self::write_itol_annotations(&path, &genome_ids, &cluster_ids, &presence, num_genomes)?;
            Some(path)
        } else {
            None
        };

        // Optionally run GrapeTree
        if Self::detect_grapetree() {
            let grapetree_prefix = downstream_dir.join("grapetree_tree");
            Self::run_grapetree(&profiles_path, &grapetree_prefix)?;
        } else {
            tracing::warn!(
                "GrapeTree is not installed. Profiles written but tree not generated. \
                 Install with: conda install -c bioconda grapetree"
            );
        }

        Ok(GrapetreeResult {
            profiles_path,
            itol_path,
            num_genomes,
            num_genes,
        })
    }

    fn name(&self) -> &str {
        "GrapeTree/iTOL Export"
    }

    fn is_available(&self) -> bool {
        // Profiles can always be exported; GrapeTree tool is optional
        true
    }

    fn required_inputs(&self) -> Vec<DownstreamInput> {
        vec![DownstreamInput::FinalGraph]
    }
}

/// Result of GrapeTree/iTOL export.
#[derive(Debug)]
pub struct GrapetreeResult {
    /// Path to the generated profiles TSV file.
    pub profiles_path: PathBuf,
    /// Path to the generated iTOL annotations file (if enabled).
    pub itol_path: Option<PathBuf>,
    /// Number of genomes in the profile.
    pub num_genomes: usize,
    /// Number of gene clusters in the profile.
    pub num_genes: usize,
}

impl DownstreamResult for GrapetreeResult {
    fn write_to(&self, _dir: &Path) -> Result<()> {
        // Files are already written by the runner at their stored paths
        if !self.profiles_path.exists() {
            return Err(Error::Output(format!(
                "Profiles file not found: {:?}",
                self.profiles_path
            )));
        }
        if let Some(ref itol_path) = self.itol_path {
            if !itol_path.exists() {
                return Err(Error::Output(format!(
                    "iTOL annotations file not found: {:?}",
                    itol_path
                )));
            }
        }
        Ok(())
    }

    fn summary(&self) -> String {
        format!(
            "GrapeTree export: {} genomes, {} gene clusters. Profiles: {}",
            self.num_genomes,
            self.num_genes,
            self.profiles_path.display()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn make_test_gml(content: &str) -> TempDir {
        let temp = TempDir::new().unwrap();
        let gml_path = temp.path().join("final_graph.gml");
        std::fs::write(&gml_path, content).unwrap();
        temp
    }

    #[test]
    fn test_parse_simple_gml() {
        let gml = r#"graph [
  directed 0
  node [
    id "c1"
    label "c1"
    support 3
    is_paralog 0
    genomes [
      "genome1"
      "genome2"
    ]
  ]
  node [
    id "c2"
    label "c2"
    support 2
    is_paralog 0
    genomes [
      "genome2"
    ]
  ]
  edge [
    source "c1"
    target "c2"
  ]
]"#;
        let temp = make_test_gml(gml);
        let graph = GrapeTreeExportRunner::parse_gml_graph(&temp.path().join("final_graph.gml")).unwrap();

        assert_eq!(graph.nodes.len(), 2);
        assert_eq!(
            graph.nodes.get(&ClusterId::new("c1")).unwrap().support,
            3
        );
        assert_eq!(
            graph.nodes.get(&ClusterId::new("c2")).unwrap().support,
            2
        );
    }

    #[test]
    fn test_build_profile_matrix() {
        let mut graph = PangenomeGraph::new();

        let mut node1 = Node::from_cluster(&{
            let mut c = crate::graph::GeneCluster::new("c1");
            c.support = 3;
            c
        });
        node1.genomes.insert(GenomeId::new("genome1"));
        node1.genomes.insert(GenomeId::new("genome2"));
        node1.genomes.insert(GenomeId::new("genome3"));
        graph.add_node(node1);

        let mut node2 = Node::from_cluster(&{
            let mut c = crate::graph::GeneCluster::new("c2");
            c.support = 2;
            c
        });
        node2.genomes.insert(GenomeId::new("genome2"));
        node2.genomes.insert(GenomeId::new("genome3"));
        graph.add_node(node2);

        let (genome_ids, cluster_ids, presence) = GrapeTreeExportRunner::build_profile_matrix(&graph);

        assert_eq!(genome_ids.len(), 3);
        assert_eq!(cluster_ids.len(), 2);
        assert_eq!(presence.len(), 3);
        assert_eq!(presence[0].len(), 2);

        // genome1 has c1 only
        let g1_idx = genome_ids.iter().position(|g| g == "genome1").unwrap();
        assert_eq!(presence[g1_idx][0], 1); // c1 present
        assert_eq!(presence[g1_idx][1], 0); // c2 absent

        // genome2 has c1 and c2
        let g2_idx = genome_ids.iter().position(|g| g == "genome2").unwrap();
        assert_eq!(presence[g2_idx][0], 1);
        assert_eq!(presence[g2_idx][1], 1);
    }

    #[test]
    fn test_profile_matrix_empty_genomes() {
        let mut graph = PangenomeGraph::new();
        let node = Node::from_cluster(&{
            let mut c = crate::graph::GeneCluster::new("c1");
            c.support = 0;
            c
        });
        graph.add_node(node);

        let (genome_ids, _cluster_ids, _presence) =
            GrapeTreeExportRunner::build_profile_matrix(&graph);

        assert!(genome_ids.is_empty());
    }

    #[test]
    fn test_grapetree_not_available() {
        // GrapeTree detection is platform/tool dependent; just verify it doesn't panic
        let available = GrapeTreeExportRunner::detect_grapetree();
        // Result depends on whether grapetree is installed on this system
        let _ = available;
    }

    #[test]
    fn test_run_grapetree_errors_when_not_installed() {
        // If GrapeTree is not installed, run_grapetree should return an error
        // (not silently return Ok as it did before the fix)
        if !GrapeTreeExportRunner::detect_grapetree() {
            let profiles = PathBuf::from("profiles.csv");
            let output = PathBuf::from("output");
            let result = GrapeTreeExportRunner::run_grapetree(&profiles, &output);
            assert!(result.is_err(), "Expected error when GrapeTree not installed, got Ok");
            let msg = result.unwrap_err().to_string();
            assert!(
                msg.to_lowercase().contains("grapetree"),
                "Error should mention grapetree, got: {}",
                msg
            );
        }
    }

    #[test]
    fn test_profile_matrix_sorted_deterministic() {
        let mut graph = PangenomeGraph::new();

        let mut node1 = Node::from_cluster(&{
            let mut c = crate::graph::GeneCluster::new("c1");
            c.support = 2;
            c
        });
        node1.genomes.insert(GenomeId::new("b_genome"));
        node1.genomes.insert(GenomeId::new("a_genome"));
        graph.add_node(node1);

        let (genome_ids, _, _) = GrapeTreeExportRunner::build_profile_matrix(&graph);

        // Should be sorted: a_genome, b_genome
        assert_eq!(genome_ids[0], "a_genome");
        assert_eq!(genome_ids[1], "b_genome");
    }
}
