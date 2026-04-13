//! Highly variable gene detection using cycle-based graph analysis.
//!
//! Implements Panaroo's `identify_possible_highly_variable` algorithm:
//! 1. Find a fundamental cycle basis via BFS spanning trees
//! 2. Filter cycles to length 3–20
//! 3. Transitively merge overlapping cycles (sharing >1 node) using union-find
//! 4. Flag nodes with support < 0.5 × max_support in merged sets with ≥5 cycles

use std::collections::{HashMap, HashSet, VecDeque};

use crate::graph::{ClusterId, PangenomeGraph};

/// Result of highly variable gene detection.
#[derive(Debug)]
pub struct HighlyVariableResult {
    /// Cluster IDs flagged as highly variable.
    pub highly_variable: HashSet<ClusterId>,
    /// Number of cycles found in the cycle basis.
    pub cycles_found: usize,
    /// Number of merged cycle sets after transitive merging.
    pub merged_sets: usize,
}

/// A merged set of overlapping cycles.
struct MergedCycleSet {
    /// All unique nodes in this merged set.
    nodes: HashSet<ClusterId>,
    /// Number of contributing cycles.
    contributing_cycles: usize,
}

/// Detect highly variable genes using Panaroo's cycle-based algorithm.
///
/// Genes in cycle-rich regions of the pangenome graph with low support
/// relative to their neighborhood are flagged as "highly variable" —
/// indicating potential hypervariable gene families (e.g., phase-variable genes,
/// antigenic variation).
pub struct HighlyVariableDetector {
    /// Minimum cycle length to consider (default: 3).
    min_cycle_length: usize,
    /// Maximum cycle length to consider (default: 20).
    max_cycle_length: usize,
    /// Minimum number of cycles in a merged set to trigger flagging (default: 5).
    min_cycles_for_flagging: usize,
    /// Support ratio threshold: nodes with support < this × max_support are flagged (default: 0.5).
    support_ratio_threshold: f64,
}

impl HighlyVariableDetector {
    /// Create a new detector with Panaroo-compatible defaults.
    pub fn new() -> Self {
        Self {
            min_cycle_length: 3,
            max_cycle_length: 20,
            min_cycles_for_flagging: 5,
            support_ratio_threshold: 0.5,
        }
    }

    /// Set minimum cycle length.
    pub fn with_min_cycle_length(mut self, n: usize) -> Self {
        self.min_cycle_length = n;
        self
    }

    /// Set maximum cycle length.
    pub fn with_max_cycle_length(mut self, n: usize) -> Self {
        self.max_cycle_length = n;
        self
    }

    /// Set minimum number of cycles for flagging.
    pub fn with_min_cycles_for_flagging(mut self, n: usize) -> Self {
        self.min_cycles_for_flagging = n;
        self
    }

    /// Set support ratio threshold.
    pub fn with_support_ratio_threshold(mut self, t: f64) -> Self {
        self.support_ratio_threshold = t;
        self
    }

    /// Detect highly variable genes in the pangenome graph.
    ///
    /// Returns the set of cluster IDs flagged as highly variable along with
    /// diagnostic statistics.
    pub fn detect(&self, graph: &PangenomeGraph) -> HighlyVariableResult {
        if graph.nodes.is_empty() {
            return HighlyVariableResult {
                highly_variable: HashSet::new(),
                cycles_found: 0,
                merged_sets: 0,
            };
        }

        // Step 1: Find cycle basis
        let cycles = self.find_cycle_basis(graph);
        let cycles_found = cycles.len();

        if cycles.is_empty() {
            return HighlyVariableResult {
                highly_variable: HashSet::new(),
                cycles_found: 0,
                merged_sets: 0,
            };
        }

        // Step 2: Merge overlapping cycles
        let merged_sets = self.merge_overlapping_cycles(&cycles);

        // Step 3: Flag highly variable genes
        let highly_variable = self.flag_highly_variable(graph, &merged_sets);

        HighlyVariableResult {
            highly_variable,
            cycles_found,
            merged_sets: merged_sets.len(),
        }
    }

    /// Find a fundamental cycle basis using BFS spanning trees.
    ///
    /// For each connected component, perform BFS to build a spanning tree.
    /// Each non-tree edge creates one fundamental cycle by tracing back
    /// through parent pointers to the lowest common ancestor.
    fn find_cycle_basis(&self, graph: &PangenomeGraph) -> Vec<HashSet<ClusterId>> {
        let mut all_cycles: Vec<HashSet<ClusterId>> = Vec::new();
        let mut global_visited: HashSet<ClusterId> = HashSet::new();

        // Process each connected component
        for start in graph.nodes.keys() {
            if global_visited.contains(start) {
                continue;
            }

            // BFS to build spanning tree
            let mut parent: HashMap<ClusterId, ClusterId> = HashMap::new();
            let mut depth: HashMap<ClusterId, usize> = HashMap::new();
            let mut bfs_order: Vec<ClusterId> = Vec::new();
            let mut queue: VecDeque<ClusterId> = VecDeque::new();

            parent.insert(start.clone(), start.clone());
            depth.insert(start.clone(), 0);
            queue.push_back(start.clone());

            while let Some(current) = queue.pop_front() {
                if global_visited.contains(&current) {
                    continue;
                }
                global_visited.insert(current.clone());
                bfs_order.push(current.clone());

                for neighbor in graph.neighbors(&current) {
                    if !parent.contains_key(neighbor) {
                        parent.insert(neighbor.clone(), current.clone());
                        depth.insert(neighbor.clone(), depth[&current] + 1);
                        queue.push_back(neighbor.clone());
                    } else if neighbor != &parent[&current] && depth.contains_key(neighbor) {
                        // Non-tree edge: trace cycle
                        // Only process each edge once (from the shallower node)
                        if depth[&current] < depth[neighbor] {
                            continue; // will be found from the deeper side
                        }

                        let cycle = self.trace_cycle(&current, neighbor, &parent, &depth);
                        let cycle_len = cycle.len();
                        if cycle_len >= self.min_cycle_length && cycle_len <= self.max_cycle_length {
                            all_cycles.push(cycle);
                        }
                    }
                }
            }
        }

        all_cycles
    }

    /// Trace a cycle from a non-tree edge through the spanning tree.
    ///
    /// Given nodes u and v connected by a non-tree edge, trace paths
    /// from u and v up through parent pointers to their lowest common
    /// ancestor, forming a cycle.
    fn trace_cycle(
        &self,
        u: &ClusterId,
        v: &ClusterId,
        parent: &HashMap<ClusterId, ClusterId>,
        depth: &HashMap<ClusterId, usize>,
    ) -> HashSet<ClusterId> {
        let mut cycle_nodes = HashSet::new();

        let mut u_curr = u.clone();
        let mut v_curr = v.clone();

        // Walk both up to equal depth
        while depth.get(&u_curr).copied().unwrap_or(0) > depth.get(&v_curr).copied().unwrap_or(0) {
            cycle_nodes.insert(u_curr.clone());
            u_curr = parent.get(&u_curr).cloned().unwrap_or(u_curr.clone());
        }
        while depth.get(&v_curr).copied().unwrap_or(0) > depth.get(&u_curr).copied().unwrap_or(0) {
            cycle_nodes.insert(v_curr.clone());
            v_curr = parent.get(&v_curr).cloned().unwrap_or(v_curr.clone());
        }

        // Walk both up until they meet at LCA
        let max_steps = depth.len() + 1; // safety bound
        let mut steps = 0;
        while u_curr != v_curr && steps < max_steps {
            cycle_nodes.insert(u_curr.clone());
            cycle_nodes.insert(v_curr.clone());
            u_curr = parent.get(&u_curr).cloned().unwrap_or(u_curr.clone());
            v_curr = parent.get(&v_curr).cloned().unwrap_or(v_curr.clone());
            steps += 1;
        }
        cycle_nodes.insert(u_curr); // LCA

        cycle_nodes
    }

    /// Merge overlapping cycles transitively using union-find.
    ///
    /// Cycles sharing more than 1 node are grouped together.
    fn merge_overlapping_cycles(&self, cycles: &[HashSet<ClusterId>]) -> Vec<MergedCycleSet> {
        let n = cycles.len();
        if n == 0 {
            return Vec::new();
        }

        // Union-Find
        let mut uf_parent: Vec<usize> = (0..n).collect();

        fn find(parent: &mut Vec<usize>, i: usize) -> usize {
            if parent[i] != i {
                parent[i] = find(parent, parent[i]);
            }
            parent[i]
        }

        fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
            let root_a = find(parent, a);
            let root_b = find(parent, b);
            if root_a != root_b {
                parent[root_b] = root_a;
            }
        }

        // Check all pairs for overlap > 1 node
        for i in 0..n {
            for j in (i + 1)..n {
                let intersection_count = cycles[i].intersection(&cycles[j]).count();
                if intersection_count > 1 {
                    union(&mut uf_parent, i, j);
                }
            }
        }

        // Group cycles by their root
        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..n {
            let root = find(&mut uf_parent, i);
            groups.entry(root).or_default().push(i);
        }

        // Build merged sets
        groups
            .into_iter()
            .map(|(_, indices)| {
                let mut all_nodes = HashSet::new();
                for &idx in &indices {
                    all_nodes.extend(cycles[idx].iter().cloned());
                }
                MergedCycleSet {
                    nodes: all_nodes,
                    contributing_cycles: indices.len(),
                }
            })
            .collect()
    }

    /// Flag highly variable genes based on merged cycle sets.
    ///
    /// For each merged set with enough contributing cycles (≥ min_cycles_for_flagging),
    /// any node whose support is below `support_ratio_threshold × max_support`
    /// in that set is flagged as highly variable.
    fn flag_highly_variable(
        &self,
        graph: &PangenomeGraph,
        merged_sets: &[MergedCycleSet],
    ) -> HashSet<ClusterId> {
        let mut highly_variable = HashSet::new();

        for set in merged_sets {
            if set.contributing_cycles < self.min_cycles_for_flagging {
                continue;
            }

            // Find max support in this merged set
            let max_support = set
                .nodes
                .iter()
                .filter_map(|id| graph.nodes.get(id))
                .map(|node| node.support)
                .max()
                .unwrap_or(0);

            if max_support == 0 {
                continue;
            }

            let threshold = (self.support_ratio_threshold * max_support as f64) as usize;

            // Flag nodes with support below threshold
            for node_id in &set.nodes {
                if let Some(node) = graph.nodes.get(node_id) {
                    if node.support < threshold {
                        highly_variable.insert(node_id.clone());
                    }
                }
            }
        }

        highly_variable
    }
}

impl Default for HighlyVariableDetector {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{GeneCluster, Node, Edge, GenomeId};

    fn make_node(id: &str, support: usize) -> Node {
        let mut cluster = GeneCluster::new(id);
        cluster.support = support;
        let mut node = Node::from_cluster(&cluster);
        node.support = support;
        node
    }

    fn make_edge(from: &str, to: &str, support: usize) -> Edge {
        Edge {
            from: ClusterId::new(from),
            to: ClusterId::new(to),
            genomes: {
                let mut g = HashSet::new();
                for i in 0..support {
                    g.insert(GenomeId::new(format!("genome_{}", i)));
                }
                g
            },
            support,
        }
    }

    #[test]
    fn test_empty_graph_no_cycles() {
        let graph = PangenomeGraph::new();
        let detector = HighlyVariableDetector::new();
        let result = detector.detect(&graph);
        assert!(result.highly_variable.is_empty());
        assert_eq!(result.cycles_found, 0);
        assert_eq!(result.merged_sets, 0);
    }

    #[test]
    fn test_simple_triangle_cycle() {
        let mut graph = PangenomeGraph::new();
        // Triangle: A -- B -- C -- A
        graph.add_node(make_node("A", 10));
        graph.add_node(make_node("B", 10));
        graph.add_node(make_node("C", 5));
        graph.add_edge(make_edge("A", "B", 10));
        graph.add_edge(make_edge("B", "C", 5));
        graph.add_edge(make_edge("C", "A", 5));

        let detector = HighlyVariableDetector::new()
            .with_min_cycles_for_flagging(1); // Lower threshold for testing
        let result = detector.detect(&graph);

        assert!(result.cycles_found >= 1, "Should find at least one cycle");
        // C has support 5, threshold = 0.5 * 10 = 5
        // 5 < 5 is false, so C should NOT be flagged
        assert!(!result.highly_variable.contains(&ClusterId::new("C")));
    }

    #[test]
    fn test_simple_triangle_with_low_support_flagged() {
        let mut graph = PangenomeGraph::new();
        // Triangle: A -- B -- C -- A where C has very low support
        graph.add_node(make_node("A", 20));
        graph.add_node(make_node("B", 20));
        graph.add_node(make_node("C", 2));
        graph.add_edge(make_edge("A", "B", 20));
        graph.add_edge(make_edge("B", "C", 2));
        graph.add_edge(make_edge("C", "A", 2));

        let detector = HighlyVariableDetector::new()
            .with_min_cycles_for_flagging(1);
        let result = detector.detect(&graph);

        // C has support 2 < 0.5 * 20 = 10, so C should be flagged
        assert!(result.highly_variable.contains(&ClusterId::new("C")));
    }

    #[test]
    fn test_linear_chain_no_cycles() {
        let mut graph = PangenomeGraph::new();
        // Linear: A -- B -- C -- D (no cycle)
        graph.add_node(make_node("A", 10));
        graph.add_node(make_node("B", 10));
        graph.add_node(make_node("C", 10));
        graph.add_node(make_node("D", 10));
        graph.add_edge(make_edge("A", "B", 10));
        graph.add_edge(make_edge("B", "C", 10));
        graph.add_edge(make_edge("C", "D", 10));

        let detector = HighlyVariableDetector::new();
        let result = detector.detect(&graph);

        assert_eq!(result.cycles_found, 0, "Linear chain should have no cycles");
        assert!(result.highly_variable.is_empty());
    }

    #[test]
    fn test_overlapping_cycles_merged() {
        let mut graph = PangenomeGraph::new();
        // Two triangles sharing edge A-B:
        //   A -- B -- C -- A (triangle 1)
        //   A -- B -- D -- A (triangle 2)
        // These share nodes A and B (2 shared nodes > 1), so should merge
        graph.add_node(make_node("A", 20));
        graph.add_node(make_node("B", 20));
        graph.add_node(make_node("C", 3));
        graph.add_node(make_node("D", 3));
        graph.add_edge(make_edge("A", "B", 20));
        graph.add_edge(make_edge("B", "C", 3));
        graph.add_edge(make_edge("C", "A", 3));
        graph.add_edge(make_edge("B", "D", 3));
        graph.add_edge(make_edge("D", "A", 3));

        let detector = HighlyVariableDetector::new()
            .with_min_cycles_for_flagging(1);
        let result = detector.detect(&graph);

        // Both C and D should be flagged (support 3 < 0.5 * 20 = 10)
        assert!(result.highly_variable.contains(&ClusterId::new("C")));
        assert!(result.highly_variable.contains(&ClusterId::new("D")));
    }

    #[test]
    fn test_non_overlapping_cycles_not_merged() {
        let mut graph = PangenomeGraph::new();
        // Two separate triangles connected by a bridge:
        //   A -- B -- C -- A (triangle 1)
        //   D -- E -- F -- D (triangle 2)
        //   C -- D (bridge, shares only 1 node with each triangle)
        graph.add_node(make_node("A", 20));
        graph.add_node(make_node("B", 20));
        graph.add_node(make_node("C", 20));
        graph.add_node(make_node("D", 20));
        graph.add_node(make_node("E", 20));
        graph.add_node(make_node("F", 20));
        graph.add_edge(make_edge("A", "B", 20));
        graph.add_edge(make_edge("B", "C", 20));
        graph.add_edge(make_edge("C", "A", 20));
        graph.add_edge(make_edge("D", "E", 20));
        graph.add_edge(make_edge("E", "F", 20));
        graph.add_edge(make_edge("F", "D", 20));
        graph.add_edge(make_edge("C", "D", 20)); // bridge

        let detector = HighlyVariableDetector::new()
            .with_min_cycles_for_flagging(5); // Need 5 cycles, only 2 found
        let result = detector.detect(&graph);

        // Not enough cycles (2 < 5), so no flags
        assert!(result.highly_variable.is_empty());
    }

    #[test]
    fn test_min_cycles_for_flagging_filter() {
        let mut graph = PangenomeGraph::new();
        // Create a graph where cycles exist but < min_cycles_for_flagging
        // Single triangle with low-support node
        graph.add_node(make_node("A", 20));
        graph.add_node(make_node("B", 20));
        graph.add_node(make_node("C", 2));
        graph.add_edge(make_edge("A", "B", 20));
        graph.add_edge(make_edge("B", "C", 2));
        graph.add_edge(make_edge("C", "A", 2));

        // Default: min_cycles_for_flagging = 5, only 1 cycle found
        let detector = HighlyVariableDetector::new();
        let result = detector.detect(&graph);
        assert!(result.highly_variable.is_empty(), "Should not flag with default threshold");
    }

    #[test]
    fn test_detector_default() {
        let detector = HighlyVariableDetector::default();
        assert_eq!(detector.min_cycle_length, 3);
        assert_eq!(detector.max_cycle_length, 20);
        assert_eq!(detector.min_cycles_for_flagging, 5);
        assert!((detector.support_ratio_threshold - 0.5).abs() < f64::EPSILON);
    }
}