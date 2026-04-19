//! Core data types for the pangenome graph.

use std::collections::{HashSet, HashMap};
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
    /// Centroid sequences (representative)
    pub centroids: Vec<Sequence>,
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
            centroids: vec![],
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
    /// Is this a highly variable gene cluster?
    pub is_highly_variable: bool,
    /// Centroid sequence (representative)
    pub centroid_sequence: Option<Sequence>,
    /// Whether this node represents a contig end
    pub is_contig_end: bool,
    /// Contig sequences where this gene appears (for missing gene recovery)
    pub contig_sequences: HashMap<String, Sequence>,
    /// Gene members per genome: genome_id -> [gene_id, gene_id, ...]
    pub gene_members: HashMap<GenomeId, Vec<String>>,
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
            is_highly_variable: false,
            centroid_sequence: cluster.centroids.first().cloned(),
            is_contig_end: false,
            contig_sequences: HashMap::new(),
            gene_members: HashMap::new(),
        }
    }

    /// Create a new node from a cluster with gene member data.
    ///
    /// Populates `gene_members` by mapping each gene ID to its genome
    /// via the `gene_data` lookup table.
    pub fn from_cluster_with_genes(
        cluster: &GeneCluster,
        gene_data: &HashMap<GeneId, Gene>,
    ) -> Self {
        let mut node = Self::from_cluster(cluster);
        for gene_id in &cluster.genes {
            if let Some(gene) = gene_data.get(gene_id) {
                node.gene_members
                    .entry(gene.genome_id.clone())
                    .or_default()
                    .push(gene_id.as_str().to_string());
            }
        }
        node
    }

    /// Add a contig sequence to this node.
    pub fn add_contig_sequence(&mut self, contig_name: impl Into<String>, sequence: Sequence) {
        self.contig_sequences.insert(contig_name.into(), sequence);
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
    /// Lookup table: gene_id -> Gene (for output writers to access contig/start/end/strand)
    pub gene_lookup: HashMap<GeneId, Gene>,
    /// Adjacency index: cluster_id -> set of neighbor cluster_ids (O(1) lookup)
    adjacency: HashMap<ClusterId, HashSet<ClusterId>>,
    /// Whether the adjacency index needs rebuilding
    adjacency_dirty: bool,
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
        // Update adjacency index incrementally
        self.adjacency.entry(edge.from.clone()).or_default().insert(edge.to.clone());
        self.adjacency.entry(edge.to.clone()).or_default().insert(edge.from.clone());
        self.edges.insert(key, edge);
    }

    /// Rebuild the adjacency index from scratch.
    pub fn rebuild_adjacency(&mut self) {
        self.adjacency.clear();
        for ((from, to), _) in &self.edges {
            self.adjacency.entry(from.clone()).or_default().insert(to.clone());
            self.adjacency.entry(to.clone()).or_default().insert(from.clone());
        }
        self.adjacency_dirty = false;
    }

    /// Get the degree of a node (number of connected edges).
    /// Uses the adjacency index for O(1) lookup.
    pub fn degree(&self, cluster_id: &ClusterId) -> usize {
        self.adjacency.get(cluster_id).map(|s| s.len()).unwrap_or(0)
    }

    /// Get neighbors of a node.
    /// Uses the adjacency index for O(1) lookup.
    pub fn neighbors(&self, cluster_id: &ClusterId) -> Vec<&ClusterId> {
        self.adjacency.get(cluster_id)
            .map(|s| s.iter().collect())
            .unwrap_or_default()
    }

    /// Check if two nodes are connected by an edge.
    pub fn has_edge(&self, from: &ClusterId, to: &ClusterId) -> bool {
        self.edges.contains_key(&(from.clone(), to.clone()))
            || self.edges.contains_key(&(to.clone(), from.clone()))
    }

    /// Remove a node and all its edges from the graph.
    pub fn remove_node(&mut self, cluster_id: &ClusterId) -> Option<Node> {
        // Remove edges and update adjacency
        if let Some(neighbors) = self.adjacency.remove(cluster_id) {
            for neighbor in &neighbors {
                if let Some(neighbor_set) = self.adjacency.get_mut(neighbor) {
                    neighbor_set.remove(cluster_id);
                }
                self.edges.remove(&(cluster_id.clone(), neighbor.clone()));
                self.edges.remove(&(neighbor.clone(), cluster_id.clone()));
            }
        }
        self.nodes.remove(cluster_id)
    }

    /// Remove an edge from the graph.
    pub fn remove_edge(&mut self, from: &ClusterId, to: &ClusterId) -> Option<Edge> {
        // Update adjacency
        if let Some(s) = self.adjacency.get_mut(from) { s.remove(to); }
        if let Some(s) = self.adjacency.get_mut(to) { s.remove(from); }
        self.edges.remove(&(from.clone(), to.clone())).or_else(|| self.edges.remove(&(to.clone(), from.clone())))
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

        let node3 = Node::from_cluster(&{
            let mut c = GeneCluster::new("c3");
            c.support = 4;
            c
        });

        graph.add_node(node1);
        graph.add_node(node2);
        graph.add_node(node3);

        let mut edge12 = Edge::new(ClusterId::new("c1"), ClusterId::new("c2"));
        edge12.add_genome(GenomeId::new("genome1"));
        graph.add_edge(edge12);

        let mut edge13 = Edge::new(ClusterId::new("c1"), ClusterId::new("c3"));
        edge13.add_genome(GenomeId::new("genome2"));
        graph.add_edge(edge13);

        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.edge_count(), 2);
        assert_eq!(graph.degree(&ClusterId::new("c1")), 2);
        assert_eq!(graph.degree(&ClusterId::new("c2")), 1);
        assert_eq!(graph.degree(&ClusterId::new("c3")), 1);
    }

    #[test]
    fn test_adjacency_index() {
        let mut graph = PangenomeGraph::new();
        graph.add_node(Node::from_cluster(&{ let mut c = GeneCluster::new("a"); c.support = 1; c }));
        graph.add_node(Node::from_cluster(&{ let mut c = GeneCluster::new("b"); c.support = 1; c }));
        graph.add_node(Node::from_cluster(&{ let mut c = GeneCluster::new("c"); c.support = 1; c }));

        graph.add_edge(Edge::new(ClusterId::new("a"), ClusterId::new("b")));
        graph.add_edge(Edge::new(ClusterId::new("a"), ClusterId::new("c")));

        assert_eq!(graph.degree(&ClusterId::new("a")), 2);
        assert_eq!(graph.degree(&ClusterId::new("b")), 1);
        assert_eq!(graph.degree(&ClusterId::new("c")), 1);

        let neighbors_a: Vec<_> = graph.neighbors(&ClusterId::new("a")).into_iter().collect();
        assert_eq!(neighbors_a.len(), 2);
    }

    #[test]
    fn test_remove_node() {
        let mut graph = PangenomeGraph::new();
        graph.add_node(Node::from_cluster(&{ let mut c = GeneCluster::new("x"); c.support = 1; c }));
        graph.add_node(Node::from_cluster(&{ let mut c = GeneCluster::new("y"); c.support = 1; c }));
        graph.add_edge(Edge::new(ClusterId::new("x"), ClusterId::new("y")));

        let removed = graph.remove_node(&ClusterId::new("x"));
        assert!(removed.is_some());
        assert_eq!(graph.node_count(), 1);
        assert_eq!(graph.edge_count(), 0);
        assert_eq!(graph.degree(&ClusterId::new("y")), 0);
    }

    #[test]
    fn test_remove_edge() {
        let mut graph = PangenomeGraph::new();
        graph.add_node(Node::from_cluster(&{ let mut c = GeneCluster::new("p"); c.support = 1; c }));
        graph.add_node(Node::from_cluster(&{ let mut c = GeneCluster::new("q"); c.support = 1; c }));
        graph.add_edge(Edge::new(ClusterId::new("p"), ClusterId::new("q")));

        let edge = graph.remove_edge(&ClusterId::new("p"), &ClusterId::new("q"));
        assert!(edge.is_some());
        assert_eq!(graph.edge_count(), 0);
        assert_eq!(graph.degree(&ClusterId::new("p")), 0);
    }

    #[test]
    fn test_node_from_cluster_with_centroid() {
        let mut cluster = GeneCluster::new("test_cluster");
        cluster.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];
        cluster.support = 3;

        let node = Node::from_cluster(&cluster);
        assert_eq!(node.centroid_sequence, Some(b"ATCGATCGATCGATCG".to_vec()));
    }

    #[test]
    fn test_gene_cluster_multiple_centroids() {
        let cluster = GeneCluster {
            id: ClusterId::new("cluster_0"),
            genes: vec![GeneId::new("gene_0")],
            centroids: vec![b"ATCG".to_vec(), b"GCTA".to_vec()],
            is_paralog: false,
            support: 1,
        };
        assert_eq!(cluster.centroids.len(), 2);
        assert_eq!(cluster.centroids[0], b"ATCG".to_vec());
        assert_eq!(cluster.centroids[1], b"GCTA".to_vec());
    }

    #[test]
    fn test_node_gene_members_default() {
        let cluster = GeneCluster::new("test_cluster");
        let node = Node::from_cluster(&cluster);
        assert!(node.gene_members.is_empty());
    }

    #[test]
    fn test_node_from_cluster_with_genes() {
        let mut cluster = GeneCluster::new("c1");
        cluster.add_gene(GeneId::new("geneA"));
        cluster.add_gene(GeneId::new("geneB"));

        let mut gene_data = std::collections::HashMap::new();
        let mut gene_a = Gene::new("geneA", GenomeId::new("genome1"));
        gene_a.contig = "contig1".to_string();
        gene_a.start = 100;
        gene_a.end = 200;
        let mut gene_b = Gene::new("geneB", GenomeId::new("genome2"));
        gene_b.contig = "contig2".to_string();

        gene_data.insert(GeneId::new("geneA"), gene_a);
        gene_data.insert(GeneId::new("geneB"), gene_b);

        let node = Node::from_cluster_with_genes(&cluster, &gene_data);
        assert_eq!(node.gene_members.len(), 2);
        assert!(node.gene_members.contains_key(&GenomeId::new("genome1")));
        assert!(node.gene_members.contains_key(&GenomeId::new("genome2")));
        assert_eq!(node.gene_members[&GenomeId::new("genome1")], vec!["geneA".to_string()]);
        assert_eq!(node.gene_members[&GenomeId::new("genome2")], vec!["geneB".to_string()]);
    }

    #[test]
    fn test_pangenome_graph_gene_lookup() {
        let mut graph = PangenomeGraph::new();
        let gene = Gene::new("geneA", GenomeId::new("genome1"));
        graph.gene_lookup.insert(GeneId::new("geneA"), gene);
        assert!(graph.gene_lookup.contains_key(&GeneId::new("geneA")));
    }
}