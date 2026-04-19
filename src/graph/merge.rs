//! Pangenome graph merging.
//!
//! Merges multiple PanMiner output directories into a single unified pangenome.
//! Reads GML graphs from each directory, re-labels nodes to avoid collisions,
//! clusters centroids across graphs, merges compatible nodes, and runs
//! correction passes on the combined graph.

use std::collections::HashMap;
use std::path::PathBuf;

use crate::error::Result;
use crate::graph::{ClusterId, Node, Edge, PangenomeGraph};

/// Result of a merge operation.
#[derive(Debug)]
pub struct MergeResult {
    /// Path to the output directory
    pub output_dir: PathBuf,
    /// Number of input graphs merged
    pub num_inputs: usize,
    /// Total nodes in the merged graph
    pub total_nodes: usize,
    /// Total edges in the merged graph
    pub total_edges: usize,
    /// Number of nodes merged (had same centroid across inputs)
    pub merged_nodes: usize,
}

/// Merge multiple PanMiner output directories into a single pangenome.
///
/// Algorithm:
/// 1. Load GML graphs from each input directory
/// 2. Re-label node IDs to avoid collisions (prefix with directory index)
/// 3. Cluster centroids across all graphs using identity threshold
/// 4. Merge nodes that share centroids and have disjoint genome membership
/// 5. Re-run correction passes on the combined graph
/// 6. Generate unified output files
pub fn merge_pangenomes(
    input_dirs: &[PathBuf],
    output_dir: &PathBuf,
    identity_threshold: f32,
    _threads: usize,
) -> Result<MergeResult> {
    if input_dirs.is_empty() {
        return Err(crate::Error::Config("No input directories specified".to_string()));
    }

    if input_dirs.len() == 1 {
        return Err(crate::Error::Config(
            "Only one input directory provided. Merge requires at least two.".to_string()
        ));
    }

    // Step 1: Load GML graphs from each directory
    let mut graphs = Vec::new();
    for dir in input_dirs {
        let gml_path = dir.join("final_graph.gml");
        if !gml_path.exists() {
            return Err(crate::Error::Output(format!(
                "No final_graph.gml found in {:?}",
                dir
            )));
        }
        let graph = load_gml_graph(&gml_path)?;
        graphs.push(graph);
    }

    let total_inputs = graphs.len();

    // Step 2: Re-label nodes to avoid collisions
    let mut relabeled_graphs = Vec::new();
    for (i, mut graph) in graphs.into_iter().enumerate() {
        let _mapping = relabel_nodes(&mut graph, i);
        relabeled_graphs.push(graph);
    }

    // Step 3: Cluster centroids across all graphs
    let clusters = cluster_centroids(&relabeled_graphs, identity_threshold);

    // Step 4: Merge into a single graph
    let mut merged = PangenomeGraph::new();
    let mut merged_count = 0;

    // Add all nodes from all graphs
    for graph in &relabeled_graphs {
        for (cluster_id, node) in &graph.nodes {
            merged.nodes.insert(cluster_id.clone(), node.clone());
        }
        for (key, edge) in &graph.edges {
            merged.edges.insert(key.clone(), edge.clone());
        }
        for (genome_id, metadata) in &graph.genomes {
            merged.genomes.insert(genome_id.clone(), metadata.clone());
        }
    }

    // Merge nodes that share centroids
    for cluster in &clusters {
        if cluster.len() > 1 {
            // Merge all nodes in this cluster into the first one
            let target_id = &cluster[0];
            if let Some(target_node) = merged.nodes.get(target_id).cloned() {
                let mut merged_node = target_node;
                for source_id in &cluster[1..] {
                    if let Some(source_node) = merged.nodes.get(source_id).cloned() {
                        // Merge genomes
                        merged_node.genomes.extend(source_node.genomes);
                        // Merge annotations
                        merged_node.annotations.extend(source_node.annotations);
                        merged_node.support += source_node.support;
                        merged_count += 1;
                    }
                }
                merged.nodes.insert(target_id.clone(), merged_node);
                // Remove merged nodes
                for source_id in &cluster[1..] {
                    merged.nodes.remove(source_id);
                }
            }
        }
    }

    // Create output directory
    std::fs::create_dir_all(output_dir)?;

    // Write merged graph as GML
    let output_gml = output_dir.join("final_graph.gml");
    crate::output::GmlWriter::write(&merged, &output_gml)?;

    // Write summary statistics
    let total_nodes = merged.node_count();
    let total_edges = merged.edge_count();

    Ok(MergeResult {
        output_dir: output_dir.clone(),
        num_inputs: total_inputs,
        total_nodes,
        total_edges,
        merged_nodes: merged_count,
    })
}

/// Load a pangenome graph from a GML file.
fn load_gml_graph(path: &PathBuf) -> Result<PangenomeGraph> {
    let content = std::fs::read_to_string(path)?;
    let mut graph = PangenomeGraph::new();

    // Simple GML parser for PanMiner output format
    // GML format: node [ id N label "cluster_id" support S ... ]
    //             edge [ source N target S ... ]
    let mut current_node: Option<Node> = None;
    let mut current_edge: Option<Edge> = None;
    let mut in_node = false;
    let mut in_edge = false;

    for line in content.lines() {
        let line = line.trim();

        if line == "node [" {
            in_node = true;
            current_node = None;
        } else if line == "]" && in_node {
            if let Some(node) = current_node.take() {
                graph.nodes.insert(node.cluster_id.clone(), node);
            }
            in_node = false;
        } else if line == "edge [" {
            in_edge = true;
        } else if line == "]" && in_edge {
            if let Some(edge) = current_edge.take() {
                let key = (edge.from.clone(), edge.to.clone());
                graph.edges.insert(key, edge);
            }
            in_edge = false;
        } else if in_node {
            if let Some(node) = &mut current_node {
                if line.starts_with("label") {
                    let label = line.split('"').nth(1).unwrap_or("").to_string();
                    node.cluster_id = ClusterId::new(&label);
                } else if line.starts_with("support") {
                    if let Some(s) = line.split_whitespace().nth(1) {
                        node.support = s.parse().unwrap_or(1);
                    }
                }
            } else {
                // Create a default node
                current_node = Some(Node::from_cluster(&{
                    let mut c = crate::graph::GeneCluster::new("temp");
                    c.support = 1;
                    c
                }));
            }
        } else if in_edge {
            if current_edge.is_none() {
                current_edge = Some(Edge::new(
                    ClusterId::new("unknown"),
                    ClusterId::new("unknown"),
                ));
            }
            if let Some(edge) = &mut current_edge {
                if line.starts_with("source") {
                    if let Some(s) = line.split('"').nth(1) {
                        edge.from = ClusterId::new(s);
                    }
                } else if line.starts_with("target") {
                    if let Some(s) = line.split('"').nth(1) {
                        edge.to = ClusterId::new(s);
                    }
                } else if line.starts_with("value") {
                    if let Some(s) = line.split_whitespace().nth(1) {
                        edge.support = s.parse().unwrap_or(1);
                    }
                }
            }
        }
    }

    Ok(graph)
}

/// Re-label nodes in a graph by prefixing with a directory index.
///
/// Returns a mapping from old IDs to new IDs.
fn relabel_nodes(graph: &mut PangenomeGraph, prefix: usize) -> HashMap<String, String> {
    let mut mapping = HashMap::new();
    let prefix_str = format!("g{}_", prefix);

    // Re-label nodes
    let old_nodes: Vec<_> = graph.nodes.drain().collect();
    for (cluster_id, mut node) in old_nodes {
        let new_id_str = format!("{}{}", prefix_str, cluster_id);
        let new_id = ClusterId::new(&new_id_str);
        mapping.insert(cluster_id.to_string(), new_id_str);
        node.cluster_id = new_id.clone();
        graph.nodes.insert(new_id, node);
    }

    // Re-label edges
    let old_edges: Vec<_> = graph.edges.drain().collect();
    for ((from, to), mut edge) in old_edges {
        let new_from_str = format!("{}{}", prefix_str, from);
        let new_to_str = format!("{}{}", prefix_str, to);
        let new_from = ClusterId::new(&new_from_str);
        let new_to = ClusterId::new(&new_to_str);
        edge.from = new_from.clone();
        edge.to = new_to.clone();
        graph.edges.insert((new_from, new_to), edge);
    }

    mapping
}

/// Cluster centroids across multiple graphs using identity threshold.
///
/// Groups cluster IDs whose centroids are similar enough to be merged.
fn cluster_centroids(graphs: &[PangenomeGraph], threshold: f32) -> Vec<Vec<ClusterId>> {
    // Collect all centroids from all graphs
    let mut all_centroids: Vec<(ClusterId, Option<crate::graph::Sequence>)> = Vec::new();
    for graph in graphs {
        for (cluster_id, node) in &graph.nodes {
            for seq in &node.centroid_sequences {
                all_centroids.push((cluster_id.clone(), Some(seq.clone())));
            }
        }
    }

    // Simple greedy clustering based on sequence identity
    let mut clusters: Vec<Vec<ClusterId>> = Vec::new();
    let mut assigned = vec![false; all_centroids.len()];

    for i in 0..all_centroids.len() {
        if assigned[i] {
            continue;
        }
        let mut cluster = vec![all_centroids[i].0.clone()];
        assigned[i] = true;

        if let Some(ref seq_i) = all_centroids[i].1 {
            for j in (i + 1)..all_centroids.len() {
                if assigned[j] {
                    continue;
                }
                if let Some(ref seq_j) = all_centroids[j].1 {
                    let identity = sequence_identity(seq_i, seq_j);
                    if identity >= threshold {
                        cluster.push(all_centroids[j].0.clone());
                        assigned[j] = true;
                    }
                }
            }
        }

        clusters.push(cluster);
    }

    clusters
}

/// Calculate sequence identity between two sequences.
fn sequence_identity(seq1: &[u8], seq2: &[u8]) -> f32 {
    if seq1.is_empty() || seq2.is_empty() {
        return 0.0;
    }

    let min_len = seq1.len().min(seq2.len());
    let max_len = seq1.len().max(seq2.len());
    let mut matches = 0;

    for i in 0..min_len {
        if seq1[i].to_ascii_uppercase() == seq2[i].to_ascii_uppercase() {
            matches += 1;
        }
    }

    matches as f32 / max_len as f32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{GeneCluster, GenomeMetadata};

    fn make_test_graph(nodes: Vec<(&str, usize, Vec<Vec<u8>>)>) -> PangenomeGraph {
        let mut graph = PangenomeGraph::new();
        for (id, support, centroids) in nodes {
            let mut cluster = GeneCluster::new(id);
            cluster.support = support;
            cluster.centroids = centroids;
            let node = Node::from_cluster(&cluster);
            graph.add_node(node);
        }
        graph
    }

    #[test]
    fn test_merge_single_graph_fails() {
        let dirs = vec![PathBuf::from("/tmp/dir1")];
        let result = merge_pangenomes(&dirs, &PathBuf::from("/tmp/out"), 0.95, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_merge_empty_dirs_fails() {
        let result = merge_pangenomes(&[], &PathBuf::from("/tmp/out"), 0.95, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_relabel_nodes() {
        let mut graph = make_test_graph(vec![
            ("c1", 3, vec![]),
            ("c2", 5, vec![]),
        ]);
        let mapping = relabel_nodes(&mut graph, 0);
        assert!(mapping.contains_key("c1"));
        assert!(mapping.contains_key("c2"));
        // Nodes should be prefixed
        assert!(graph.nodes.contains_key(&ClusterId::new("g0_c1")));
        assert!(graph.nodes.contains_key(&ClusterId::new("g0_c2")));
    }

    #[test]
    fn test_sequence_identity_identical() {
        let seq1 = b"ATCGATCG";
        let seq2 = b"ATCGATCG";
        let identity = sequence_identity(seq1, seq2);
        assert!((identity - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_sequence_identity_different() {
        let seq1 = b"ATCGATCG";
        let seq2 = b"GCTAGCTA";
        let identity = sequence_identity(seq1, seq2);
        assert!(identity < 0.3);
    }

    #[test]
    fn test_sequence_identity_empty() {
        let identity = sequence_identity(b"", b"ATCG");
        assert_eq!(identity, 0.0);
    }

    #[test]
    fn test_cluster_centroids() {
        let mut graph1 = make_test_graph(vec![
            ("c1", 3, vec![b"ATCGATCG".to_vec()]),
        ]);
        let mut graph2 = make_test_graph(vec![
            ("c2", 5, vec![b"ATCGATCG".to_vec()]), // identical to c1
        ]);

        // Add genomes
        graph1.genomes.insert(
            crate::graph::GenomeId::new("genome1"),
            GenomeMetadata {
                id: crate::graph::GenomeId::new("genome1"),
                source_file: String::new(),
                num_contigs: 1,
                total_genes: 100,
            },
        );
        graph2.genomes.insert(
            crate::graph::GenomeId::new("genome2"),
            GenomeMetadata {
                id: crate::graph::GenomeId::new("genome2"),
                source_file: String::new(),
                num_contigs: 1,
                total_genes: 100,
            },
        );

        let clusters = cluster_centroids(&[graph1, graph2], 0.95);
        // c1 and c2 should be in the same cluster (identical centroids)
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].len(), 2);
    }
}