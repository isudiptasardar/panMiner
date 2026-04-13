//! Gene neighborhood extraction from pangenome graphs.
//!
//! This module provides the [`GeneNeighborhoodExtractor`] which extracts genes
//! within N hops of a seed gene in the pangenome graph using BFS traversal.
//!
//! # Algorithm
//!
//! 1. Parse GML from `output_dir / "final_graph.gml"` -> `PangenomeGraph`
//! 2. Find the seed cluster node by `cluster_id`
//! 3. BFS from seed node to `max_depth` hops, recording hop distance per reachable node
//! 4. Collect subgraph edges among visited nodes
//! 5. Write `neighborhood_genes.csv` (cluster_id, support, annotation, hop_distance, num_genomes)
//! 6. Write `neighborhood_subgraph.gml` (Cytoscape-compatible, with hop_distance as node attribute)
//!
//! # Example
//!
//! ```no_run
//! use std::path::Path;
//! use panminer::downstream::exploration::GeneNeighborhoodExtractor;
//! use panminer::downstream::DownstreamResult;
//! use panminer::graph::ClusterId;
//!
//! # fn main() -> panminer::Result<()> {
//! let extractor = GeneNeighborhoodExtractor::with_default_depth(ClusterId::new("cluster_001"));
//! let result = extractor.run(Path::new("panminer_output"))?;
//! result.write_to(Path::new("neighborhood_results"))?;
//! # Ok(())
//! # }
//! ```

use std::collections::{HashMap, HashSet, VecDeque};
use std::io::Write;
use std::path::Path;

use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};
use crate::error::Result;
use crate::graph::{ClusterId, Edge, Node, PangenomeGraph};

/// GeneNeighborhoodExtractor extracts genes within N hops of a seed gene.
pub struct GeneNeighborhoodExtractor {
    /// Cluster ID to start BFS from (seed gene)
    seed_gene: ClusterId,
    /// Maximum BFS depth (default: 5)
    max_depth: usize,
}

impl GeneNeighborhoodExtractor {
    /// Create a new neighborhood extractor.
    ///
    /// # Arguments
    ///
    /// * `seed_gene` - Cluster ID of the seed gene to start BFS from
    /// * `max_depth` - Maximum number of hops to traverse (default: 5)
    pub fn new(seed_gene: ClusterId, max_depth: usize) -> Self {
        Self {
            seed_gene,
            max_depth,
        }
    }

    /// Create with default max depth of 5.
    pub fn with_default_depth(seed_gene: ClusterId) -> Self {
        Self {
            seed_gene,
            max_depth: 5,
        }
    }

    /// Extract neighborhood from a GML file in the output directory.
    pub fn run(&self, output_dir: &Path) -> Result<NeighborhoodResult> {
        let gml_path = output_dir.join("final_graph.gml");
        let graph = parse_gml_graph(&gml_path)?;

        let distances = bfs_neighborhood(&graph, &self.seed_gene, self.max_depth);
        let visited: HashSet<ClusterId> = distances.keys().cloned().collect();

        let mut visited_nodes = Vec::new();
        for (cluster_id, &hop_distance) in &distances {
            if let Some(node) = graph.nodes.get(cluster_id) {
                visited_nodes.push(NeighborhoodNode {
                    cluster_id: cluster_id.to_string(),
                    support: node.support,
                    annotation: node.annotations.iter().next().cloned(),
                    hop_distance,
                    num_genomes: node.genomes.len(),
                    is_paralog: node.is_paralog,
                    is_highly_variable: node.is_highly_variable,
                });
            }
        }

        let edges = collect_subgraph_edges(&graph, &visited);

        Ok(NeighborhoodResult {
            seed_gene: self.seed_gene.to_string(),
            max_depth: self.max_depth,
            visited_nodes,
            edges,
        })
    }
}

impl DownstreamRunner for GeneNeighborhoodExtractor {
    type Output = NeighborhoodResult;

    fn run(&self, output_dir: &Path) -> Result<Self::Output> {
        self.run(output_dir)
    }

    fn name(&self) -> &str {
        "GeneNeighborhoodExtractor"
    }

    fn is_available(&self) -> bool {
        true // Native Rust, always available
    }

    fn required_inputs(&self) -> Vec<DownstreamInput> {
        vec![DownstreamInput::FinalGraph]
    }
}

/// A node in the neighborhood with its hop distance from the seed.
#[derive(Debug, Clone)]
pub struct NeighborhoodNode {
    /// Cluster ID
    pub cluster_id: String,
    /// Number of genomes containing this cluster
    pub support: usize,
    /// Gene annotation (e.g., product name)
    pub annotation: Option<String>,
    /// Hop distance from the seed gene
    pub hop_distance: usize,
    /// Number of genomes in the cluster
    pub num_genomes: usize,
    /// Whether this is a paralog cluster
    pub is_paralog: bool,
    /// Whether this is a highly variable gene cluster
    pub is_highly_variable: bool,
}

/// Result of neighborhood extraction.
pub struct NeighborhoodResult {
    /// The seed gene cluster ID
    pub seed_gene: String,
    /// Maximum BFS depth used
    pub max_depth: usize,
    /// All visited nodes with their hop distances
    pub visited_nodes: Vec<NeighborhoodNode>,
    /// Edges among the visited nodes (from, to)
    pub edges: Vec<(String, String)>,
}

impl DownstreamResult for NeighborhoodResult {
    fn write_to(&self, dir: &Path) -> Result<()> {
        std::fs::create_dir_all(dir)?;

        // Write neighborhood_genes.csv
        let csv_path = dir.join("neighborhood_genes.csv");
        let mut wtr = csv::Writer::from_path(&csv_path)?;
        wtr.write_record(&["cluster_id", "support", "annotation", "hop_distance", "num_genomes", "is_paralog", "is_highly_variable"])?;
        for node in &self.visited_nodes {
            let paralog_str = if node.is_paralog { "true" } else { "false" };
            let hv_str = if node.is_highly_variable { "true" } else { "false" };
            wtr.write_record(&[
                node.cluster_id.as_str(),
                &node.support.to_string(),
                node.annotation.as_deref().unwrap_or(""),
                &node.hop_distance.to_string(),
                &node.num_genomes.to_string(),
                paralog_str,
                hv_str,
            ])?;
        }
        wtr.flush()?;

        // Write neighborhood_subgraph.gml
        let gml_path = dir.join("neighborhood_subgraph.gml");
        write_neighborhood_gml(&self.visited_nodes, &self.edges, &gml_path)?;

        Ok(())
    }

    fn summary(&self) -> String {
        format!(
            "GeneNeighborhoodExtractor: seed={}, max_depth={}, visited_nodes={}, edges={}",
            self.seed_gene,
            self.max_depth,
            self.visited_nodes.len(),
            self.edges.len()
        )
    }
}

impl DownstreamResult for GeneNeighborhoodExtractor {
    fn write_to(&self, _dir: &Path) -> Result<()> {
        // No-op: use run() to get NeighborhoodResult
        Ok(())
    }

    fn summary(&self) -> String {
        format!(
            "GeneNeighborhoodExtractor: seed={}, max_depth={}",
            self.seed_gene, self.max_depth
        )
    }
}

/// Parse GML file into a PangenomeGraph.
///
/// This follows the simple line-by-line pattern from `src/graph/merge.rs`.
fn parse_gml_graph(path: &Path) -> Result<PangenomeGraph> {
    let content = std::fs::read_to_string(path)?;
    let mut graph = PangenomeGraph::new();

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
            current_edge = None;
        } else if line == "]" && in_edge {
            if let Some(edge) = current_edge.take() {
                let key = (edge.from.clone(), edge.to.clone());
                graph.edges.insert(key, edge);
            }
            in_edge = false;
        } else if in_node {
            if current_node.is_none() {
                current_node = Some(Node::from_cluster(&{
                    let mut c = crate::graph::GeneCluster::new("temp");
                    c.support = 1;
                    c
                }));
            }
            if let Some(node) = &mut current_node {
                if line.starts_with("label") {
                    let label = line.split('"').nth(1).unwrap_or("").to_string();
                    node.cluster_id = ClusterId::new(&label);
                } else if line.starts_with("support") {
                    if let Some(s) = line.split_whitespace().nth(1) {
                        node.support = s.parse().unwrap_or(1);
                    }
                } else if line.starts_with("is_paralog") {
                    if let Some(s) = line.split_whitespace().nth(1) {
                        node.is_paralog = s == "1";
                    }
                } else if line.starts_with("is_highly_variable") {
                    if let Some(s) = line.split_whitespace().nth(1) {
                        node.is_highly_variable = s == "1";
                    }
                } else if line.starts_with("annotation") {
                    if let Some(ann) = line.split('"').nth(1) {
                        node.annotations.insert(ann.to_string());
                    }
                } else if line.starts_with("genomes") {
                    // genomes are listed after the node, parse them
                    // Format: genomes number_value or genomes [ ... ]
                    if let Some(s) = line.split_whitespace().nth(1) {
                        let count: usize = s.parse().unwrap_or(0);
                        // Add placeholder genomes based on support
                        for i in 0..count {
                            node.genomes
                                .insert(crate::graph::GenomeId::new(format!("genome_{}", i)));
                        }
                    }
                }
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
                } else if line.starts_with("genomes") {
                    if let Some(s) = line.split_whitespace().nth(1) {
                        let count: usize = s.parse().unwrap_or(0);
                        for i in 0..count {
                            edge.genomes
                                .insert(crate::graph::GenomeId::new(format!("genome_{}", i)));
                        }
                        edge.support = edge.genomes.len();
                    }
                }
            }
        }
    }

    Ok(graph)
}

/// Perform BFS from seed node, returning map of cluster_id -> hop_distance.
fn bfs_neighborhood(
    graph: &PangenomeGraph,
    seed: &ClusterId,
    max_depth: usize,
) -> HashMap<ClusterId, usize> {
    let mut distances: HashMap<ClusterId, usize> = HashMap::new();
    let mut queue: VecDeque<(ClusterId, usize)> = VecDeque::new();

    // Initialize with seed at depth 0
    if graph.nodes.contains_key(seed) {
        distances.insert(seed.clone(), 0);
        queue.push_back((seed.clone(), 0));
    }

    while let Some((current, depth)) = queue.pop_front() {
        if depth >= max_depth {
            continue;
        }

        // Get neighbors
        for neighbor in graph.neighbors(&current) {
            if !distances.contains_key(neighbor) {
                distances.insert(neighbor.clone(), depth + 1);
                queue.push_back((neighbor.clone(), depth + 1));
            }
        }
    }

    distances
}

/// Collect edges among visited nodes (both endpoints in visited set).
fn collect_subgraph_edges(
    graph: &PangenomeGraph,
    visited: &HashSet<ClusterId>,
) -> Vec<(String, String)> {
    let mut edges = Vec::new();
    for ((from, to), _) in &graph.edges {
        if visited.contains(from) && visited.contains(to) {
            edges.push((from.to_string(), to.to_string()));
        }
    }
    edges
}

/// Write neighborhood subgraph as GML with hop_distance node attribute.
fn write_neighborhood_gml(
    nodes: &[NeighborhoodNode],
    edges: &[(String, String)],
    path: &Path,
) -> Result<()> {
    let mut file = std::fs::File::create(path)?;

    writeln!(&mut file, "graph [")?;
    writeln!(&mut file, "  directed 0")?;

    // Write nodes with hop_distance attribute
    for node in nodes {
        writeln!(&mut file, "  node [")?;
        writeln!(&mut file, "    id \"{}\"", node.cluster_id)?;
        writeln!(&mut file, "    label \"{}\"", node.cluster_id)?;
        writeln!(&mut file, "    support {}", node.support)?;
        writeln!(&mut file, "    is_paralog {}", if node.is_paralog { 1 } else { 0 })?;
        writeln!(&mut file, "    is_highly_variable {}", if node.is_highly_variable { 1 } else { 0 })?;
        writeln!(&mut file, "    hop_distance {}", node.hop_distance)?;
        writeln!(&mut file, "    num_genomes {}", node.num_genomes)?;
        if let Some(ref ann) = node.annotation {
            let escaped = ann.replace('"', "\\\"");
            writeln!(&mut file, "    annotation \"{}\"", escaped)?;
        }
        writeln!(&mut file, "  ]")?;
    }

    // Write edges
    for (from, to) in edges {
        writeln!(&mut file, "  edge [")?;
        writeln!(&mut file, "    source \"{}\"", from)?;
        writeln!(&mut file, "    target \"{}\"", to)?;
        writeln!(&mut file, "  ]")?;
    }

    writeln!(&mut file, "]")?;

    Ok(())
}

/// Escape a string for GML format.
#[allow(dead_code)]
fn escape_gml(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{GeneCluster, GenomeId};

    fn make_test_graph() -> PangenomeGraph {
        let mut graph = PangenomeGraph::new();

        // Create nodes: c1 -> c2 -> c3 -> c4 (linear chain)
        // c1 also connects to c5
        for (id, support) in [("c1", 3), ("c2", 5), ("c3", 4), ("c4", 2), ("c5", 6)] {
            let node = Node::from_cluster(&{
                let mut c = GeneCluster::new(id);
                c.support = support;
                c
            });
            graph.add_node(node);
        }

        // Add edges
        for (from, to) in [("c1", "c2"), ("c2", "c3"), ("c3", "c4"), ("c1", "c5")] {
            let mut edge = Edge::new(ClusterId::new(from), ClusterId::new(to));
            edge.add_genome(GenomeId::new("g1"));
            graph.add_edge(edge);
        }

        graph
    }

    #[test]
    fn test_bfs_depth_0() {
        let graph = make_test_graph();
        let distances = bfs_neighborhood(&graph, &ClusterId::new("c1"), 0);
        // Depth 0 should only include the seed
        assert!(distances.contains_key(&ClusterId::new("c1")));
        assert_eq!(distances.len(), 1);
    }

    #[test]
    fn test_bfs_depth_1() {
        let graph = make_test_graph();
        let distances = bfs_neighborhood(&graph, &ClusterId::new("c1"), 1);
        // c1 -> c2, c1 -> c5 at depth 1
        assert_eq!(distances.get(&ClusterId::new("c1")), Some(&0));
        assert_eq!(distances.get(&ClusterId::new("c2")), Some(&1));
        assert_eq!(distances.get(&ClusterId::new("c5")), Some(&1));
        assert!(!distances.contains_key(&ClusterId::new("c3")));
    }

    #[test]
    fn test_bfs_depth_2() {
        let graph = make_test_graph();
        let distances = bfs_neighborhood(&graph, &ClusterId::new("c1"), 2);
        // c1 -> c2 -> c3 at depth 2
        assert_eq!(distances.get(&ClusterId::new("c3")), Some(&2));
        assert!(!distances.contains_key(&ClusterId::new("c4")));
    }

    #[test]
    fn test_bfs_depth_3() {
        let graph = make_test_graph();
        let distances = bfs_neighborhood(&graph, &ClusterId::new("c1"), 3);
        // c1 -> c2 -> c3 -> c4 at depth 3
        assert_eq!(distances.get(&ClusterId::new("c4")), Some(&3));
        // All nodes should be reachable
        assert_eq!(distances.len(), 5);
    }

    #[test]
    fn test_collect_subgraph_edges() {
        let graph = make_test_graph();
        let visited: HashSet<_> = [
            ClusterId::new("c1"),
            ClusterId::new("c2"),
            ClusterId::new("c3"),
        ].into_iter().collect();

        let edges = collect_subgraph_edges(&graph, &visited);
        // Only c1-c2 and c2-c3 have both endpoints in visited
        assert_eq!(edges.len(), 2);
        assert!(edges.contains(&("c1".to_string(), "c2".to_string())));
        assert!(edges.contains(&("c2".to_string(), "c3".to_string())));
    }

    #[test]
    fn test_neighborhood_result_write() {
        let result = NeighborhoodResult {
            seed_gene: "c1".to_string(),
            max_depth: 3,
            visited_nodes: vec![
                NeighborhoodNode {
                    cluster_id: "c1".to_string(),
                    support: 3,
                    annotation: Some("hypothetical protein".to_string()),
                    hop_distance: 0,
                    num_genomes: 2,
                    is_paralog: false,
                    is_highly_variable: false,
                },
            ],
            edges: vec![("c1".to_string(), "c2".to_string())],
        };

        let dir = std::env::temp_dir().join("neighborhood_test");
        result.write_to(&dir).unwrap();

        let csv_path = dir.join("neighborhood_genes.csv");
        assert!(csv_path.exists());

        let gml_path = dir.join("neighborhood_subgraph.gml");
        assert!(gml_path.exists());

        // Cleanup
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_seed_not_in_graph() {
        let graph = make_test_graph();
        let distances = bfs_neighborhood(&graph, &ClusterId::new("nonexistent"), 3);
        // No nodes reachable from nonexistent seed
        assert!(distances.is_empty());
    }
}