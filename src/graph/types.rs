//! Core data types for the pangenome graph.

use std::collections::HashSet;
use std::hash::Hash;
use serde::{Serialize, Deserialize};

/// Gene sequence (nucleotide bytes).
pub type Sequence = Vec<u8>;

/// Unique identifier for a gene.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct GeneId(String);

impl GeneId {
    pub fn new(id: impl Into<String>) -> Self {
        Self(id.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::fmt::Display for GeneId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Unique identifier for a genome.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct GenomeId(String);

impl GenomeId {
    pub fn new(id: impl Into<String>) -> Self {
        Self(id.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::fmt::Display for GenomeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Unique identifier for a gene cluster.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct ClusterId(String);

impl ClusterId {
    pub fn new(id: impl Into<String>) -> Self {
        Self(id.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::fmt::Display for ClusterId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Strand orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Default for Strand {
    fn default() -> Self {
        Self::Unknown
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

/// A gene from a genome.
#[derive(Debug, Clone)]
pub struct Gene {
    /// Unique gene identifier
    pub id: GeneId,
    /// Gene sequence (nucleotides)
    pub sequence: Sequence,
    /// Genome this gene belongs to
    pub genome_id: GenomeId,
    /// Contig/chromosome name
    pub contig: String,
    /// Start position (1-based)
    pub start: usize,
    /// End position (1-based, inclusive)
    pub end: usize,
    /// Strand orientation
    pub strand: Strand,
    /// Gene annotation (e.g., product name)
    pub annotation: Option<String>,
}

impl Gene {
    /// Create a new gene with the given ID.
    pub fn new(id: impl Into<String>, genome_id: GenomeId) -> Self {
        Self {
            id: GeneId::new(id),
            sequence: Vec::new(),
            genome_id,
            contig: String::new(),
            start: 0,
            end: 0,
            strand: Strand::Unknown,
            annotation: None,
        }
    }

    /// Get the length of the gene.
    pub fn length(&self) -> usize {
        if self.end >= self.start {
            self.end - self.start + 1
        } else {
            0
        }
    }
}

/// A cluster of orthologous genes (COG).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneCluster {
    /// Unique cluster identifier
    pub id: ClusterId,
    /// Genes in this cluster
    pub genes: Vec<GeneId>,
    /// Centroid sequence (representative)
    pub centroid: Option<Sequence>,
    /// Whether this cluster contains paralogs
    pub is_paralog: bool,
    /// Number of genomes containing this cluster
    pub support: usize,
}

impl GeneCluster {
    /// Create a new empty cluster.
    pub fn new(id: impl Into<String>) -> Self {
        Self {
            id: ClusterId::new(id),
            genes: Vec::new(),
            centroid: None,
            is_paralog: false,
            support: 0,
        }
    }

    /// Add a gene to this cluster.
    pub fn add_gene(&mut self, gene_id: GeneId) {
        self.genes.push(gene_id);
    }

    /// Check if this cluster contains a specific gene.
    pub fn contains(&self, gene_id: &GeneId) -> bool {
        self.genes.contains(gene_id)
    }

    /// Get the number of genes in this cluster.
    pub fn len(&self) -> usize {
        self.genes.len()
    }

    /// Check if this cluster is empty.
    pub fn is_empty(&self) -> bool {
        self.genes.is_empty()
    }
}

/// Edge key for graph edges.
pub type EdgeKey = (ClusterId, ClusterId);

/// A node in the pangenome graph.
#[derive(Debug, Clone)]
pub struct Node {
    /// Cluster ID for this node
    pub cluster_id: ClusterId,
    /// Number of genomes containing this cluster
    pub support: usize,
    /// Genomes that have this cluster
    pub genomes: HashSet<GenomeId>,
    /// Gene annotations in this cluster
    pub annotations: HashSet<String>,
    /// Is this a paralog cluster?
    pub is_paralog: bool,
    /// Centroid sequence (representative)
    pub centroid_sequence: Option<Sequence>,
}

impl Node {
    /// Create a new node from a cluster.
    pub fn from_cluster(cluster: &GeneCluster) -> Self {
        Self {
            cluster_id: cluster.id.clone(),
            support: cluster.support,
            genomes: HashSet::new(),
            annotations: HashSet::new(),
            is_paralog: cluster.is_paralog,
            centroid_sequence: cluster.centroid.clone(),
        }
    }
}

/// An edge in the pangenome graph.
#[derive(Debug, Clone)]
pub struct Edge {
    /// Source cluster ID
    pub from: ClusterId,
    /// Target cluster ID
    pub to: ClusterId,
    /// Genomes that have this adjacency
    pub genomes: HashSet<GenomeId>,
    /// Support count (number of genomes)
    pub support: usize,
}

impl Edge {
    /// Create a new edge.
    pub fn new(from: ClusterId, to: ClusterId) -> Self {
        Self {
            from,
            to,
            genomes: HashSet::new(),
            support: 0,
        }
    }

    /// Add a genome to this edge.
    pub fn add_genome(&mut self, genome_id: GenomeId) {
        self.genomes.insert(genome_id);
        self.support = self.genomes.len();
    }
}

/// Metadata about a genome.
#[derive(Debug, Clone)]
pub struct GenomeMetadata {
    /// Genome ID
    pub id: GenomeId,
    /// Source file path
    pub source_file: String,
    /// Number of contigs
    pub num_contigs: usize,
    /// Total gene count
    pub total_genes: usize,
}

/// The pangenome graph structure.
#[derive(Debug, Clone, Default)]
pub struct PangenomeGraph {
    /// Nodes indexed by cluster ID
    pub nodes: std::collections::HashMap<ClusterId, Node>,
    /// Edges indexed by (from, to) tuple
    pub edges: std::collections::HashMap<EdgeKey, Edge>,
    /// Metadata for each genome
    pub genomes: std::collections::HashMap<GenomeId, GenomeMetadata>,
}

impl PangenomeGraph {
    /// Create a new empty pangenome graph.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a node to the graph.
    pub fn add_node(&mut self, node: Node) {
        self.nodes.insert(node.cluster_id.clone(), node);
    }

    /// Add an edge to the graph.
    pub fn add_edge(&mut self, edge: Edge) {
        let key = (edge.from.clone(), edge.to.clone());
        self.edges.insert(key, edge);
    }

    /// Get the degree of a node (number of connected edges).
    pub fn degree(&self, cluster_id: &ClusterId) -> usize {
        self.edges
            .iter()
            .filter(|((from, to), _)| from == cluster_id || to == cluster_id)
            .count()
    }

    /// Get neighbors of a node.
    pub fn neighbors(&self, cluster_id: &ClusterId) -> Vec<&ClusterId> {
        self.edges
            .iter()
            .filter_map(|((from, to), _)| {
                if from == cluster_id {
                    Some(to)
                } else if to == cluster_id {
                    Some(from)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Get the number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Get the number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Check if the graph is empty.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gene_id() {
        let id = GeneId::new("gene_001");
        assert_eq!(id.as_str(), "gene_001");
        assert_eq!(id.to_string(), "gene_001");
    }

    #[test]
    fn test_gene_cluster() {
        let mut cluster = GeneCluster::new("cluster_001");
        cluster.add_gene(GeneId::new("gene_001"));
        cluster.add_gene(GeneId::new("gene_002"));

        assert_eq!(cluster.len(), 2);
        assert!(cluster.contains(&GeneId::new("gene_001")));
        assert!(!cluster.is_empty());
    }

    #[test]
    fn test_pangenome_graph() {
        let mut graph = PangenomeGraph::new();

        let node1 = Node::from_cluster(&{
            let mut c = GeneCluster::new("c1");
            c.support = 5;
            c
        });

        let node2 = Node::from_cluster(&{
            let mut c = GeneCluster::new("c2");
            c.support = 3;
            c
        });

        graph.add_node(node1);
        graph.add_node(node2);

        let mut edge = Edge::new(ClusterId::new("c1"), ClusterId::new("c2"));
        edge.add_genome(GenomeId::new("genome1"));
        graph.add_edge(edge);

        assert_eq!(graph.node_count(), 2);
        assert_eq!(graph.edge_count(), 1);
        assert_eq!(graph.degree(&ClusterId::new("c1")), 1);
    }

    #[test]
    fn test_node_from_cluster_with_centroid() {
        let mut cluster = GeneCluster::new("test_cluster");
        cluster.centroid = Some(b"ATCGATCGATCGATCG".to_vec());
        cluster.support = 3;

        let node = Node::from_cluster(&cluster);
        assert_eq!(node.centroid_sequence, Some(b"ATCGATCGATCGATCG".to_vec()));
    }
}