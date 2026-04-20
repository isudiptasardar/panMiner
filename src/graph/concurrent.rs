//! Concurrent graph structure using DashMap for thread-safe updates.

use dashmap::DashMap;
use rayon::prelude::*;
use std::collections::HashSet;
use std::sync::Mutex;

use super::types::{ClusterId, Edge, EdgeKey, GenomeId, Node, PangenomeGraph};
use crate::io::PartialGraph;

/// Thread-safe pangenome graph using DashMap.
///
/// Concurrent reads and single-key writes (add_node, add_edge_genome) are
/// safe via DashMap's internal sharded locking. However, multi-step mutations
/// (merge_nodes, remove_node) that touch multiple DashMap entries must be
/// serialized through `merge_lock` to prevent race conditions.
pub struct ConcurrentGraph {
    /// Nodes indexed by cluster ID
    pub nodes: DashMap<ClusterId, Node>,
    /// Edges indexed by (from, to) tuple
    pub edges: DashMap<EdgeKey, Edge>,
    /// Adjacency index: cluster_id -> set of neighbor cluster_ids (O(1) lookup)
    adjacency: DashMap<ClusterId, HashSet<ClusterId>>,
    /// Serializes multi-step mutations (merge, remove) to prevent race conditions
    merge_lock: Mutex<()>,
}

impl ConcurrentGraph {
    /// Create a new empty concurrent graph.
    pub fn new() -> Self {
        Self {
            nodes: DashMap::new(),
            edges: DashMap::new(),
            adjacency: DashMap::new(),
            merge_lock: Mutex::new(()),
        }
    }

    /// Create a concurrent graph with pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            nodes: DashMap::with_capacity(capacity),
            edges: DashMap::with_capacity(capacity * 2), // Usually more edges than nodes
            adjacency: DashMap::with_capacity(capacity),
            merge_lock: Mutex::new(()),
        }
    }

    /// Add a node to the graph (thread-safe).
    pub fn add_node(&self, node: Node) {
        self.nodes.insert(node.cluster_id.clone(), node);
    }

    /// Add an edge to the graph (thread-safe, undirected).
    pub fn add_edge(&self, edge: Edge) {
        let key = if edge.from < edge.to {
            (edge.from.clone(), edge.to.clone())
        } else {
            (edge.to.clone(), edge.from.clone())
        };
        // Update adjacency index incrementally
        self.adjacency.entry(key.0.clone()).or_insert_with(HashSet::new).insert(key.1.clone());
        self.adjacency.entry(key.1.clone()).or_insert_with(HashSet::new).insert(key.0.clone());
        self.edges.insert(key, edge);
    }

    /// Add or update an edge with a genome (thread-safe, undirected).
    ///
    /// If the edge already exists, the genome is added to the edge's genome set.
    pub fn add_edge_genome(&self, from: ClusterId, to: ClusterId, genome: GenomeId) {
        let key = if from < to {
            (from.clone(), to.clone())
        } else {
            (to.clone(), from.clone())
        };

        // Update adjacency index incrementally
        self.adjacency.entry(key.0.clone()).or_insert_with(HashSet::new).insert(key.1.clone());
        self.adjacency.entry(key.1.clone()).or_insert_with(HashSet::new).insert(key.0.clone());

        self.edges
            .entry(key)
            .and_modify(|existing| {
                existing.genomes.insert(genome.clone());
                existing.support = existing.genomes.len();
            })
            .or_insert_with(|| {
                // Ensure the Edge struct also stores the nodes in consistent order
                let mut edge = if from < to {
                    Edge::new(from, to)
                } else {
                    Edge::new(to, from)
                };
                edge.add_genome(genome);
                edge
            });
    }

    /// Get the degree of a node (number of connected edges).
    /// Uses the adjacency index for O(1) lookup.
    pub fn degree(&self, cluster_id: &ClusterId) -> usize {
        self.adjacency.get(cluster_id).map(|s| s.len()).unwrap_or(0)
    }

    /// Get the neighbors of a node.
    /// Uses the adjacency index for O(1) lookup.
    pub fn neighbors(&self, cluster_id: &ClusterId) -> Vec<ClusterId> {
        self.adjacency.get(cluster_id)
            .map(|s| s.value().iter().cloned().collect())
            .unwrap_or_default()
    }

    /// Check if a node has degree 1 (only one connected edge).
    pub fn is_degree_one(&self, cluster_id: &ClusterId) -> bool {
        self.degree(cluster_id) == 1
    }

    /// Find all low-support nodes (for contamination removal).
    ///
    /// A node is considered low-support if:
    /// - It has support below the threshold
    /// - AND it has degree <= 1 (only one connected edge)
    pub fn find_low_support_nodes(&self, threshold: usize) -> Vec<ClusterId> {
        self.nodes
            .iter()
            .par_bridge()
            .filter(|entry| {
                let node = entry.value();
                node.support < threshold && self.is_degree_one(entry.key())
            })
            .map(|entry| entry.key().clone())
            .collect()
    }

    /// Remove a node and its connected edges (thread-safe).
    /// Uses the adjacency index for O(degree) instead of O(E).
    /// Acquires merge_lock to serialize multi-step mutations.
    pub fn remove_node(&self, cluster_id: &ClusterId) {
        let _guard = self.merge_lock.lock().unwrap();

        // Get neighbors from adjacency (O(1) lookup)
        let neighbors: Vec<ClusterId> = {
            if let Some(adj) = self.adjacency.get(cluster_id) {
                adj.value().iter().cloned().collect()
            } else {
                Vec::new()
            }
        };

        // Remove connected edges (O(degree) instead of O(E))
        for neighbor in &neighbors {
            let key = if cluster_id < neighbor {
                (cluster_id.clone(), neighbor.clone())
            } else {
                (neighbor.clone(), cluster_id.clone())
            };
            self.edges.remove(&key);
        }

        // Remove this node from all neighbors' adjacency sets
        for neighbor in &neighbors {
            if let Some(mut adj) = self.adjacency.get_mut(neighbor) {
                adj.remove(cluster_id);
            }
        }

        // Remove node's own adjacency entry
        self.adjacency.remove(cluster_id);

        // Remove node
        self.nodes.remove(cluster_id);
    }

    /// Remove an edge from the graph (thread-safe).
    /// Updates the adjacency index for both endpoints.
    /// Acquires merge_lock to serialize multi-step mutations.
    pub fn remove_edge(&self, from: &ClusterId, to: &ClusterId) -> Option<Edge> {
        let _guard = self.merge_lock.lock().unwrap();

        // Update adjacency: remove each endpoint from the other's set
        if let Some(mut adj) = self.adjacency.get_mut(from) {
            adj.remove(to);
        }
        if let Some(mut adj) = self.adjacency.get_mut(to) {
            adj.remove(from);
        }

        // Remove edge from edges DashMap (canonical key ordering: min, max)
        let key = if from < to {
            (from.clone(), to.clone())
        } else {
            (to.clone(), from.clone())
        };
        self.edges.remove(&key).map(|(_, edge)| edge)
    }

    /// Remove multiple nodes sequentially under the merge lock.
    ///
    /// This must be sequential (not parallel) because each removal
    /// is a multi-step operation that must be serialized to prevent
    /// race conditions on the adjacency index.
    pub fn remove_nodes_parallel(&self, nodes: &[ClusterId]) {
        // Each remove_node acquires merge_lock internally
        for node in nodes {
            self.remove_node(node);
        }
    }

    /// Merge the source node into the target node.
    ///
    /// This rewires all edges connected to `source` to connect to `target` instead,
    /// merges the node data (support and annotations), and then deletes `source`.
    /// This preserves the contiguity of the graph during error correction.
    /// Uses the adjacency index for O(degree) instead of O(E).
    /// Acquires merge_lock to prevent race conditions on multi-step mutations.
    pub fn merge_nodes(&self, target: &ClusterId, source: &ClusterId) {
        if target == source {
            return;
        }

        let _guard = self.merge_lock.lock().unwrap();

        // Get source's neighbors from adjacency (O(1))
        let source_neighbors: Vec<ClusterId> = {
            if let Some(adj) = self.adjacency.get(source) {
                adj.value().iter().cloned().collect()
            } else {
                Vec::new()
            }
        };

        // Remove all edges connected to source and collect genomes for rewiring (O(degree))
        let mut edges_to_rewire = Vec::new();
        for neighbor in &source_neighbors {
            let key = if source < neighbor {
                (source.clone(), neighbor.clone())
            } else {
                (neighbor.clone(), source.clone())
            };

            if let Some((_, edge)) = self.edges.remove(&key) {
                if neighbor != target {
                    edges_to_rewire.push((neighbor.clone(), edge.genomes.clone()));
                }
                // If neighbor == target, the edge is just deleted (not rewired)
            }
        }

        // Update adjacency: remove source from all neighbors' adjacency sets
        for neighbor in &source_neighbors {
            if let Some(mut adj) = self.adjacency.get_mut(neighbor) {
                adj.remove(source);
            }
        }

        // Remove source's own adjacency entry
        self.adjacency.remove(source);

        // Rewire edges to target (add_edge_genome updates adjacency)
        for (neighbor, genomes) in edges_to_rewire {
            for genome in genomes {
                self.add_edge_genome(target.clone(), neighbor.clone(), genome);
            }
        }

        // Merge node data
        if let Some((_, source_node)) = self.nodes.remove(source) {
            self.nodes.entry(target.clone()).and_modify(|target_node| {
                target_node.support += source_node.support;
                target_node.annotations.extend(source_node.annotations);
                target_node.is_paralog |= source_node.is_paralog;
                target_node.is_highly_variable |= source_node.is_highly_variable;
                target_node.centroid_sequences.extend(source_node.centroid_sequences);
                target_node.contig_end_genomes.extend(source_node.contig_end_genomes);
                // Merge gene members from source into target
                for (genome_id, gene_ids) in source_node.gene_members {
                    target_node.gene_members.entry(genome_id).or_default().extend(gene_ids);
                }
                // Merge genomes set
                target_node.genomes.extend(source_node.genomes);
            });
        }
    }

    /// Merge multiple partial graphs into this graph.
    pub fn merge_from(&self, partials: Vec<PartialGraph>) {
        // Merge edges from all partial graphs in parallel
        partials
            .par_iter()
            .flat_map(|partial| {
                // Convert adjacencies to edges
                // Tuple format: (genome_id, from_cluster, to_cluster)
                partial.adjacencies.par_iter().map(|(genome_id, from, to)| {
                    (ClusterId::new(from), ClusterId::new(to), GenomeId::new(genome_id))
                })
            })
            .for_each(|(from, to, genome)| {
                self.add_edge_genome(from, to, genome);
            });
    }

    /// Get the number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Get the number of edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Convert to a standard PangenomeGraph.
    pub fn to_standard(&self) -> PangenomeGraph {
        let mut graph = PangenomeGraph::new();

        for entry in self.nodes.iter() {
            graph.add_node(entry.value().clone());
        }

        for entry in self.edges.iter() {
            graph.add_edge(entry.value().clone());
        }

        graph
    }
}

impl Default for ConcurrentGraph {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{GeneCluster, ClusterId, GenomeId, Node};

    #[test]
    fn test_concurrent_graph_basic() {
        let graph = ConcurrentGraph::new();

        let node = Node::from_cluster(&{
            let mut c = GeneCluster::new("c1");
            c.support = 5;
            c
        });

        graph.add_node(node);
        assert_eq!(graph.node_count(), 1);
    }

    #[test]
    fn test_concurrent_edge_update() {
        let graph = ConcurrentGraph::new();

        // Add same edge multiple times with different genomes
        graph.add_edge_genome(ClusterId::new("c1"), ClusterId::new("c2"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("c1"), ClusterId::new("c2"), GenomeId::new("g2"));
        graph.add_edge_genome(ClusterId::new("c1"), ClusterId::new("c2"), GenomeId::new("g3"));

        // Should be one edge with 3 genomes
        assert_eq!(graph.edge_count(), 1);
        let edge = graph.edges.get(&(ClusterId::new("c1"), ClusterId::new("c2"))).unwrap();
        assert_eq!(edge.support, 3);
    }

    #[test]
    fn test_find_low_support_nodes() {
        let graph = ConcurrentGraph::with_capacity(100);

        // Add a low-support node
        let low_node = Node::from_cluster(&{
            let mut c = GeneCluster::new("low");
            c.support = 1;
            c
        });
        graph.add_node(low_node);

        // Add a high-support node
        let high_node = Node::from_cluster(&{
            let mut c = GeneCluster::new("high");
            c.support = 100;
            c
        });
        graph.add_node(high_node);

        // Add edge between them (makes low node degree-1)
        graph.add_edge_genome(ClusterId::new("low"), ClusterId::new("high"), GenomeId::new("g1"));

        let low_support = graph.find_low_support_nodes(2);
        assert_eq!(low_support.len(), 1);
        assert_eq!(low_support[0], ClusterId::new("low"));
    }

    #[test]
    fn test_undirected_edges() {
        let graph = ConcurrentGraph::new();

        let c1 = ClusterId::new("c1");
        let c2 = ClusterId::new("c2");

        // Add edge c1 -> c2
        graph.add_edge_genome(c1.clone(), c2.clone(), GenomeId::new("g1"));

        // Add edge c2 -> c1 (reverse contig)
        graph.add_edge_genome(c2.clone(), c1.clone(), GenomeId::new("g2"));

        // Should be treated as a single undirected edge
        assert_eq!(graph.edge_count(), 1);

        // Support should be 2
        let key = if c1 < c2 { (c1, c2) } else { (c2, c1) };
        let edge = graph.edges.get(&key).unwrap();
        assert_eq!(edge.support, 2);
    }

    #[test]
    fn test_merge_nodes() {
        let graph = ConcurrentGraph::new();

        // Path: x -> a -> b -> y
        graph.add_node(Node::from_cluster(&{
            let mut c = GeneCluster::new("x"); c.support = 1; c
        }));
        graph.add_node(Node::from_cluster(&{
            let mut c = GeneCluster::new("a"); c.support = 1; c
        }));
        graph.add_node(Node::from_cluster(&{
            let mut c = GeneCluster::new("b"); c.support = 1; c
        }));
        graph.add_node(Node::from_cluster(&{
            let mut c = GeneCluster::new("y"); c.support = 1; c
        }));

        let g1 = GenomeId::new("g1");
        graph.add_edge_genome(ClusterId::new("x"), ClusterId::new("a"), g1.clone());
        graph.add_edge_genome(ClusterId::new("a"), ClusterId::new("b"), g1.clone());
        graph.add_edge_genome(ClusterId::new("b"), ClusterId::new("y"), g1.clone());

        // Merge b into a
        graph.merge_nodes(&ClusterId::new("a"), &ClusterId::new("b"));

        // b should be gone
        assert!(graph.nodes.get(&ClusterId::new("b")).is_none());

        // y should now be connected to a
        let key_a_y = if ClusterId::new("a") < ClusterId::new("y") {
            (ClusterId::new("a"), ClusterId::new("y"))
        } else {
            (ClusterId::new("y"), ClusterId::new("a"))
        };

        assert!(graph.edges.get(&key_a_y).is_some());
        assert_eq!(graph.degree(&ClusterId::new("a")), 2); // connected to x and y
    }
}