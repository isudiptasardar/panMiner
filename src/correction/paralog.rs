//! Paralog resolution with shortest-path distance and synteny context.
//!
//! Detects paralog clusters (same genome appears multiple times) and resolves
//! them using shortest-path distance as the primary method (matching Panaroo's
//! `nx.shortest_path_length`), with BFS context vector similarity (depth 5)
//! as fallback when no path exists within `max_context` depth.
//!
//! Paralogs are flagged during clustering via `is_paralog`, then resolved here
//! by comparing their graph neighborhoods. Paralog copies with similar synteny
//! context are merged together.

use std::collections::{HashMap, HashSet, VecDeque};

use crate::error::Result;
use crate::graph::ConcurrentGraph;
use crate::graph::{ClusterId, GenomeId, GeneCluster};

/// Compute shortest path distance between two nodes using BFS.
///
/// Returns `Some(distance)` if a path exists within `max_depth` hops,
/// or `None` if no path is found. Matching Panaroo's `nx.shortest_path_length`
/// as the primary paralog resolution method.
fn shortest_path_distance(
    graph: &ConcurrentGraph,
    from: &ClusterId,
    to: &ClusterId,
    max_depth: usize,
) -> Option<usize> {
    if from == to {
        return Some(0);
    }
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();
    visited.insert(from.clone());
    queue.push_back((from.clone(), 0usize));

    while let Some((current, depth)) = queue.pop_front() {
        if depth >= max_depth {
            continue;
        }
        for neighbor in graph.neighbors(&current) {
            if neighbor == *to {
                return Some(depth + 1);
            }
            if visited.insert(neighbor.clone()) {
                queue.push_back((neighbor, depth + 1));
            }
        }
    }
    None
}

/// Resolves paralog clusters by comparing synteny context
/// and merging copies that share similar neighborhoods.
pub struct ParalogResolver {
    /// BFS depth for context vector computation (default: 5, matching Panaroo)
    max_context: usize,
}

impl ParalogResolver {
    /// Create a new paralog resolver with default settings.
    pub fn new() -> Self {
        Self { max_context: 5 }
    }

    /// Set the BFS depth for context vector computation.
    pub fn with_max_context(mut self, depth: usize) -> Self {
        self.max_context = depth;
        self
    }

    /// Detect and resolve paralog clusters in the graph.
    ///
    /// Finds all nodes marked as paralog and merges copies that share
    /// similar synteny context (neighboring gene structure).
    ///
    /// This must run BEFORE mistranslation correction (Panaroo Step 4),
    /// because paralogs need distinct nodes so that neighbor-based
    /// comparison doesn't confuse paralog copies.
    pub fn resolve(&self, graph: &ConcurrentGraph) -> Result<ParalogStats> {
        // Step 1: Find paralog nodes
        let paralog_ids: Vec<ClusterId> = graph
            .nodes
            .iter()
            .filter(|e| e.value().is_paralog)
            .map(|e| e.key().clone())
            .collect();

        let paralogs_detected = paralog_ids.len();

        if paralog_ids.is_empty() {
            tracing::info!("Paralog resolution: no paralog clusters found");
            return Ok(ParalogStats {
                paralogs_detected: 0,
                nodes_merged: 0,
            });
        }

        tracing::info!(
            "Paralog resolution: found {} paralog nodes",
            paralogs_detected
        );

        // Step 2: Compute context vectors for all paralog nodes
        let context_vectors: HashMap<ClusterId, ContextVector> = paralog_ids
            .iter()
            .map(|id| (id.clone(), self.compute_context_vector(graph, id)))
            .collect();

        // Step 3: Merge paralog copies that share similar synteny context
        let merged_count = self.merge_by_context(graph, &paralog_ids, &context_vectors)?;

        tracing::info!(
            "Paralog resolution: {} paralogs detected, {} merged",
            paralogs_detected,
            merged_count
        );

        Ok(ParalogStats {
            paralogs_detected,
            nodes_merged: merged_count,
        })
    }

    /// Merge paralog copies that share the same synteny context.
    ///
    /// Uses shortest-path distance as the primary resolution method (matching
    /// Panaroo's `nx.shortest_path_length`), falling back to BFS context vector
    /// similarity when no path exists within `max_context` depth. Nodes with
    /// similar neighborhoods (score >= 0.5) and similar-length centroid sequences
    /// (within 20%) are merged together.
    fn merge_by_context(
        &self,
        graph: &ConcurrentGraph,
        paralog_ids: &[ClusterId],
        context_vectors: &HashMap<ClusterId, ContextVector>,
    ) -> Result<usize> {
        let mut to_merge: Vec<(ClusterId, ClusterId)> = Vec::new();

        for i in 0..paralog_ids.len() {
            for j in (i + 1)..paralog_ids.len() {
                let id_a = &paralog_ids[i];
                let id_b = &paralog_ids[j];

                // Skip if either was already removed by a prior merge
                if graph.nodes.get(id_a).is_none() || graph.nodes.get(id_b).is_none() {
                    continue;
                }

                let node_a = graph.nodes.get(id_a).unwrap();
                let node_b = graph.nodes.get(id_b).unwrap();

                // Quick length check: only merge paralogs with similar-length sequences
                let seq_a = node_a.centroid_sequences.first().map(|s| s.as_slice()).unwrap_or(&[]);
                let seq_b = node_b.centroid_sequences.first().map(|s| s.as_slice()).unwrap_or(&[]);
                let len_a = seq_a.len().max(1);
                let len_b = seq_b.len().max(1);
                let len_ratio = len_a.min(len_b) as f64 / len_a.max(len_b) as f64;
                if len_ratio < 0.8 {
                    continue;
                }

                // Primary: shortest-path distance (Panaroo's nx.shortest_path_length)
                // Fallback: context vector similarity when no path exists
                let similarity = if let Some(distance) =
                    shortest_path_distance(graph, id_a, id_b, self.max_context)
                {
                    1.0 / (1.0 + distance as f64)
                } else if let (Some(ctx_a), Some(ctx_b)) =
                    (context_vectors.get(id_a), context_vectors.get(id_b))
                {
                    ctx_a.similarity(ctx_b)
                } else {
                    continue;
                };

                // Merge if similarity exceeds threshold
                if similarity >= 0.5 {
                    // Merge lower-support node into higher-support node
                    if node_a.support >= node_b.support {
                        to_merge.push((id_a.clone(), id_b.clone()));
                    } else {
                        to_merge.push((id_b.clone(), id_a.clone()));
                    }
                }
            }
        }

        // Execute merges
        let mut merged_count = 0;
        for (target, source) in &to_merge {
            if graph.nodes.get(target).is_some() && graph.nodes.get(source).is_some() {
                graph.merge_nodes(target, source);
                merged_count += 1;
            }
        }

        Ok(merged_count)
    }

    /// Compute a context vector for a node using BFS to the given depth.
    ///
    /// The context vector captures the neighborhood structure around a node,
    /// represented as a map from neighbor cluster IDs to their BFS distance.
    /// Uses the adjacency index for O(degree) neighbor lookups.
    fn compute_context_vector(&self, graph: &ConcurrentGraph, start: &ClusterId) -> ContextVector {
        let mut visited = HashMap::new();
        let mut queue = VecDeque::new();

        visited.insert(start.clone(), 0);
        queue.push_back(start.clone());

        while let Some(current) = queue.pop_front() {
            let depth = visited.get(&current).copied().unwrap_or(0);
            if depth >= self.max_context {
                continue;
            }

            for neighbor in graph.neighbors(&current) {
                if !visited.contains_key(&neighbor) {
                    visited.insert(neighbor.clone(), depth + 1);
                    queue.push_back(neighbor);
                }
            }
        }

        // Remove the start node itself from the context
        visited.remove(start);

        ContextVector { neighbors: visited }
    }

    /// Mark clusters that contain paralogs during graph construction.
    ///
    /// A cluster is a paralog if the same genome appears more than once
    /// (i.e., two genes from the same genome were assigned to the same cluster).
    pub fn mark_paralog_clusters(clusters: &mut [GeneCluster], genes: &[crate::graph::Gene]) {
        let gene_to_genome: HashMap<String, GenomeId> = genes
            .iter()
            .map(|g| (g.id.to_string(), g.genome_id.clone()))
            .collect();

        for cluster in clusters.iter_mut() {
            let mut seen_genomes = HashSet::new();
            let mut has_paralog = false;

            for gene_id in &cluster.genes {
                if let Some(genome) = gene_to_genome.get(gene_id.as_str()) {
                    if seen_genomes.contains(genome) {
                        has_paralog = true;
                        break;
                    }
                    seen_genomes.insert(genome.clone());
                }
            }

            cluster.is_paralog = has_paralog;
        }
    }
}

impl Default for ParalogResolver {
    fn default() -> Self {
        Self::new()
    }
}

/// Context vector for synteny comparison.
///
/// Maps neighboring cluster IDs to their BFS distance from the start node.
/// Used to compare whether two paralog nodes share similar neighborhoods,
/// following Panaroo's approach of using context vectors for paralog resolution.
#[derive(Debug, Clone)]
struct ContextVector {
    neighbors: HashMap<ClusterId, usize>,
}

impl ContextVector {
    /// Compute similarity between two context vectors.
    ///
    /// Uses Panaroo's formula: sum(1 / (1 + |depth_a - depth_b|)) for
    /// shared neighbors, normalized by the total number of unique neighbors.
    fn similarity(&self, other: &ContextVector) -> f64 {
        let all_keys: HashSet<&ClusterId> = self.neighbors.keys().chain(other.neighbors.keys()).collect();

        if all_keys.is_empty() {
            return 0.0;
        }

        let mut score = 0.0;
        for key in &all_keys {
            let depth_a = self.neighbors.get(*key);
            let depth_b = other.neighbors.get(*key);

            match (depth_a, depth_b) {
                (Some(a), Some(b)) => {
                    // Both have this neighbor — similarity inversely proportional
                    // to the difference in BFS depth
                    score += 1.0 / (1.0 + (*a as f64 - *b as f64).abs());
                }
                (Some(_), None) | (None, Some(_)) => {
                    // Only one has this neighbor — penalty (0 contribution)
                    // but still counts towards normalization
                }
                _ => {}
            }
        }

        // Normalize by total unique neighbors
        score / all_keys.len() as f64
    }
}

/// Statistics from paralog resolution.
#[derive(Debug, Clone)]
pub struct ParalogStats {
    /// Number of paralog nodes detected
    pub paralogs_detected: usize,
    /// Number of paralog groups merged by context
    pub nodes_merged: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Node, Gene, GeneId, Strand};

    #[test]
    fn test_paralog_resolver_creation() {
        let resolver = ParalogResolver::new();
        assert_eq!(resolver.max_context, 5);
    }

    #[test]
    fn test_context_vector_similarity() {
        // Two nodes with identical neighborhoods
        let ctx_a = ContextVector {
            neighbors: [
                (ClusterId::new("n1"), 1),
                (ClusterId::new("n2"), 2),
                (ClusterId::new("n3"), 3),
            ]
            .iter()
            .cloned()
            .collect(),
        };

        let ctx_b = ContextVector {
            neighbors: [
                (ClusterId::new("n1"), 1),
                (ClusterId::new("n2"), 2),
                (ClusterId::new("n3"), 3),
            ]
            .iter()
            .cloned()
            .collect(),
        };

        // Identical context → similarity = 1.0
        let sim = ctx_a.similarity(&ctx_b);
        assert!(sim > 0.9, "Expected similarity > 0.9, got {}", sim);
    }

    #[test]
    fn test_context_vector_no_overlap() {
        // Two nodes with completely different neighborhoods
        let ctx_a = ContextVector {
            neighbors: [
                (ClusterId::new("a1"), 1),
                (ClusterId::new("a2"), 2),
            ]
            .iter()
            .cloned()
            .collect(),
        };

        let ctx_b = ContextVector {
            neighbors: [
                (ClusterId::new("b1"), 1),
                (ClusterId::new("b2"), 2),
            ]
            .iter()
            .cloned()
            .collect(),
        };

        // No overlap → similarity = 0
        let sim = ctx_a.similarity(&ctx_b);
        assert!(sim < 0.1, "Expected similarity < 0.1, got {}", sim);
    }

    #[test]
    fn test_context_vector_partial_overlap() {
        // Two nodes with partial overlap
        let ctx_a = ContextVector {
            neighbors: [
                (ClusterId::new("n1"), 1),
                (ClusterId::new("n2"), 2),
            ]
            .iter()
            .cloned()
            .collect(),
        };

        let ctx_b = ContextVector {
            neighbors: [
                (ClusterId::new("n1"), 1),  // shared, same depth
                (ClusterId::new("n3"), 2),   // different
            ]
            .iter()
            .cloned()
            .collect(),
        };

        // 1 shared neighbor out of 3 total, depth match = 1.0
        // similarity = 1.0 / 3 ≈ 0.33
        let sim = ctx_a.similarity(&ctx_b);
        assert!(
            sim > 0.2 && sim < 0.5,
            "Expected similarity between 0.2 and 0.5, got {}",
            sim
        );
    }

    #[test]
    fn test_mark_paralog_clusters() {
        // Cluster with two genes from the same genome → paralog
        let mut cluster = GeneCluster::new("c1");
        cluster.add_gene(GeneId::new("g1"));
        cluster.add_gene(GeneId::new("g2"));

        let genes = vec![
            Gene {
                id: GeneId::new("g1"),
                sequence: Vec::new(),
                genome_id: GenomeId::new("genome1"),
                contig: "contig1".to_string(),
                start: 1,
                end: 100,
                strand: Strand::Plus,
                annotation: None,
            },
            Gene {
                id: GeneId::new("g2"),
                sequence: Vec::new(),
                genome_id: GenomeId::new("genome1"), // Same genome → paralog
                contig: "contig1".to_string(),
                start: 200,
                end: 300,
                strand: Strand::Plus,
                annotation: None,
            },
        ];

        let mut clusters = vec![cluster];
        ParalogResolver::mark_paralog_clusters(&mut clusters, &genes);
        assert!(clusters[0].is_paralog, "Cluster should be marked as paralog");
    }

    #[test]
    fn test_mark_non_paralog_clusters() {
        // Cluster with genes from different genomes → not paralog
        let mut cluster = GeneCluster::new("c1");
        cluster.add_gene(GeneId::new("g1"));
        cluster.add_gene(GeneId::new("g2"));

        let genes = vec![
            Gene {
                id: GeneId::new("g1"),
                sequence: Vec::new(),
                genome_id: GenomeId::new("genome1"),
                contig: "contig1".to_string(),
                start: 1,
                end: 100,
                strand: Strand::Plus,
                annotation: None,
            },
            Gene {
                id: GeneId::new("g2"),
                sequence: Vec::new(),
                genome_id: GenomeId::new("genome2"), // Different genome
                contig: "contig1".to_string(),
                start: 1,
                end: 100,
                strand: Strand::Plus,
                annotation: None,
            },
        ];

        let mut clusters = vec![cluster];
        ParalogResolver::mark_paralog_clusters(&mut clusters, &genes);
        assert!(!clusters[0].is_paralog, "Cluster should NOT be marked as paralog");
    }

    #[test]
    fn test_compute_context_vector() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Build: start -> a -> b -> c
        let nodes = vec![
            Node::from_cluster(&{ let mut c = GeneCluster::new("start"); c.support = 5; c }),
            Node::from_cluster(&{ let mut c = GeneCluster::new("a"); c.support = 3; c }),
            Node::from_cluster(&{ let mut c = GeneCluster::new("b"); c.support = 3; c }),
            Node::from_cluster(&{ let mut c = GeneCluster::new("c"); c.support = 3; c }),
        ];
        for node in nodes {
            graph.add_node(node);
        }

        graph.add_edge_genome(ClusterId::new("start"), ClusterId::new("a"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("a"), ClusterId::new("b"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("b"), ClusterId::new("c"), GenomeId::new("g1"));

        let resolver = ParalogResolver::new();
        let ctx = resolver.compute_context_vector(&graph, &ClusterId::new("start"));

        // "start" should have neighbors: a(1), b(2), c(3)
        assert_eq!(ctx.neighbors.len(), 3);
        assert_eq!(ctx.neighbors.get(&ClusterId::new("a")), Some(&1));
        assert_eq!(ctx.neighbors.get(&ClusterId::new("b")), Some(&2));
        assert_eq!(ctx.neighbors.get(&ClusterId::new("c")), Some(&3));
    }

    #[test]
    fn test_resolve_no_paralogs() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Non-paralog node
        let node = Node::from_cluster(&{ let mut c = GeneCluster::new("a"); c.support = 5; c });
        graph.add_node(node);

        let resolver = ParalogResolver::new();
        let stats = resolver.resolve(&graph).unwrap();

        assert_eq!(stats.paralogs_detected, 0);
        assert_eq!(stats.nodes_merged, 0);
    }

    #[test]
    fn test_resolve_merges_similar_paralogs() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Two paralog nodes with similar neighborhoods should be merged.
        // With shortest-path as primary: direct edge gives distance 1,
        // similarity = 1/(1+1) = 0.5, meeting the merge threshold.
        let mut para1 = GeneCluster::new("para1");
        para1.support = 3;
        para1.is_paralog = true;
        para1.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];

        let mut para2 = GeneCluster::new("para2");
        para2.support = 2;
        para2.is_paralog = true;
        para2.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];

        // Shared neighbors: x and y
        let mut node_x = GeneCluster::new("x");
        node_x.support = 5;
        let mut node_y = GeneCluster::new("y");
        node_y.support = 5;

        graph.add_node(Node::from_cluster(&para1));
        graph.add_node(Node::from_cluster(&para2));
        graph.add_node(Node::from_cluster(&node_x));
        graph.add_node(Node::from_cluster(&node_y));

        // Both paralogs connected to x and y (similar context)
        graph.add_edge_genome(ClusterId::new("para1"), ClusterId::new("x"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("para1"), ClusterId::new("y"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("para2"), ClusterId::new("x"), GenomeId::new("g2"));
        graph.add_edge_genome(ClusterId::new("para2"), ClusterId::new("y"), GenomeId::new("g2"));
        // Direct edge between paralogs: shortest path = 1, similarity = 0.5
        graph.add_edge_genome(ClusterId::new("para1"), ClusterId::new("para2"), GenomeId::new("g3"));

        let resolver = ParalogResolver::new();
        let stats = resolver.resolve(&graph).unwrap();

        assert_eq!(stats.paralogs_detected, 2);
        // Should merge because direct neighbors (distance 1, similarity 0.5)
        assert_eq!(stats.nodes_merged, 1);
        // After merge, should have 3 nodes (para1 absorbed para2)
        assert_eq!(graph.node_count(), 3);
    }

    #[test]
    fn test_shortest_path_distance() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Build linear graph: A -- B -- C
        for id in ["A", "B", "C"] {
            graph.add_node(Node::from_cluster(&{
                let mut c = GeneCluster::new(id);
                c.support = 1;
                c
            }));
        }
        graph.add_edge_genome(ClusterId::new("A"), ClusterId::new("B"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("B"), ClusterId::new("C"), GenomeId::new("g1"));

        // A -> C is 2 hops
        assert_eq!(
            shortest_path_distance(&graph, &ClusterId::new("A"), &ClusterId::new("C"), 5),
            Some(2)
        );
        // A -> B is 1 hop
        assert_eq!(
            shortest_path_distance(&graph, &ClusterId::new("A"), &ClusterId::new("B"), 5),
            Some(1)
        );
        // A -> A is 0 (same node)
        assert_eq!(
            shortest_path_distance(&graph, &ClusterId::new("A"), &ClusterId::new("A"), 5),
            Some(0)
        );
        // D does not exist in the graph, so no path
        assert_eq!(
            shortest_path_distance(&graph, &ClusterId::new("A"), &ClusterId::new("D"), 5),
            None
        );
        // C -> A is also 2 (undirected)
        assert_eq!(
            shortest_path_distance(&graph, &ClusterId::new("C"), &ClusterId::new("A"), 5),
            Some(2)
        );
    }

    #[test]
    fn test_shortest_path_distance_depth_limit() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Build: A -- B -- C -- D
        for id in ["A", "B", "C", "D"] {
            graph.add_node(Node::from_cluster(&{
                let mut c = GeneCluster::new(id);
                c.support = 1;
                c
            }));
        }
        graph.add_edge_genome(ClusterId::new("A"), ClusterId::new("B"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("B"), ClusterId::new("C"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("C"), ClusterId::new("D"), GenomeId::new("g1"));

        // A -> D is 3 hops, reachable with max_depth=5
        assert_eq!(
            shortest_path_distance(&graph, &ClusterId::new("A"), &ClusterId::new("D"), 5),
            Some(3)
        );
        // Not reachable with max_depth=2 (BFS stops before reaching D)
        assert_eq!(
            shortest_path_distance(&graph, &ClusterId::new("A"), &ClusterId::new("D"), 2),
            None
        );
    }

    #[test]
    fn test_shortest_path_resolves_paralogs() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Build graph where para1 and para2 are connected through a short path:
        //   para1 -- bridge -- para2
        // Shortest path distance = 2, so similarity = 1/(1+2) ≈ 0.33 < 0.5
        // They should NOT be merged via shortest path alone.
        //
        // But if they share a direct neighbor (distance = 1):
        //   para1 -- shared -- para2
        // Then similarity = 1/(1+1) = 0.5, meeting the threshold.
        let mut para1 = GeneCluster::new("para1");
        para1.support = 3;
        para1.is_paralog = true;
        para1.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];

        let mut para2 = GeneCluster::new("para2");
        para2.support = 2;
        para2.is_paralog = true;
        para2.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];

        let mut shared = GeneCluster::new("shared");
        shared.support = 5;

        graph.add_node(Node::from_cluster(&para1));
        graph.add_node(Node::from_cluster(&para2));
        graph.add_node(Node::from_cluster(&shared));

        // Both paralogs directly connected to shared (distance = 1 each via shared)
        // para1 -> shared -> para2: distance = 2, similarity = 1/3 ≈ 0.33
        graph.add_edge_genome(ClusterId::new("para1"), ClusterId::new("shared"), GenomeId::new("g1"));
        graph.add_edge_genome(ClusterId::new("para2"), ClusterId::new("shared"), GenomeId::new("g2"));

        let resolver = ParalogResolver::new();
        let stats = resolver.resolve(&graph).unwrap();

        assert_eq!(stats.paralogs_detected, 2);
        // Shortest path distance = 2, similarity = 0.33 < 0.5, so no merge
        // Context vectors also will not help here (only 1 shared neighbor, 1 unique each)
        assert_eq!(stats.nodes_merged, 0);
    }

    #[test]
    fn test_shortest_path_direct_neighbor_paralogs() {
        let graph = ConcurrentGraph::with_capacity(10);

        // Build graph where para1 and para2 are directly connected (distance = 1):
        //   para1 -- para2
        // similarity = 1/(1+1) = 0.5, meeting the threshold.
        let mut para1 = GeneCluster::new("para1");
        para1.support = 3;
        para1.is_paralog = true;
        para1.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];

        let mut para2 = GeneCluster::new("para2");
        para2.support = 2;
        para2.is_paralog = true;
        para2.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];

        graph.add_node(Node::from_cluster(&para1));
        graph.add_node(Node::from_cluster(&para2));

        // Direct edge: shortest path = 1, similarity = 0.5
        graph.add_edge_genome(ClusterId::new("para1"), ClusterId::new("para2"), GenomeId::new("g1"));

        let resolver = ParalogResolver::new();
        let stats = resolver.resolve(&graph).unwrap();

        assert_eq!(stats.paralogs_detected, 2);
        // Direct neighbors: distance 1, similarity 0.5 >= 0.5, should merge
        assert_eq!(stats.nodes_merged, 1);
        assert_eq!(graph.node_count(), 1);
    }
}