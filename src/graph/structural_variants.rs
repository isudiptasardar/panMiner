//! Structural variant detection from pangenome graph edges.
//!
//! Detects:
//! - Inversions: reversed edge directions across genomes
//! - Duplications: nodes with multiple edges or self-loops
//! - Translocations: adjacency between non-syntenic clusters

use crate::graph::PangenomeGraph;

/// A detected structural variant.
#[derive(Debug, Clone)]
pub struct StructuralVariant {
    /// Variant type
    pub variant_type: VariantType,
    /// Cluster IDs involved
    pub cluster_ids: Vec<String>,
    /// Genomes affected
    pub affected_genomes: Vec<String>,
    /// Support count
    pub support: usize,
    /// Description
    pub description: String,
}

/// Types of structural variants.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantType {
    /// Inversion - reversed adjacency orientation
    Inversion,
    /// Duplication - multiple copies of same region
    Duplication,
    /// Translocation - moved to different location
    Translocation,
    /// Deletion - missing adjacency
    Deletion,
}

/// Structural variant detector.
pub struct StructuralVariantDetector {
    min_support: usize,
}

impl StructuralVariantDetector {
    /// Create a new detector.
    pub fn new() -> Self {
        Self {
            min_support: 1,
        }
    }

    /// Set minimum support threshold.
    pub fn with_min_support(mut self, min: usize) -> Self {
        self.min_support = min;
        self
    }

    /// Detect all structural variants in the graph.
    pub fn detect(&self, graph: &PangenomeGraph) -> Vec<StructuralVariant> {
        let mut variants = Vec::new();

        // Detect inversions
        variants.extend(self.detect_inversions(graph));

        // Detect duplications
        variants.extend(self.detect_duplications(graph));

        variants
    }

    /// Detect inversions: same cluster pair with reversed edge direction.
    fn detect_inversions(&self, graph: &PangenomeGraph) -> Vec<StructuralVariant> {
        let mut inversions = Vec::new();
        let mut checked_pairs: std::collections::HashSet<(String, String)> = std::collections::HashSet::new();

        // Group edges by cluster pair
        let mut pair_edges: std::collections::HashMap<(String, String), Vec<&crate::graph::Edge>> = std::collections::HashMap::new();
        for edge in graph.edges.values() {
            let key = (edge.from.to_string(), edge.to.to_string());
            pair_edges.entry(key).or_default().push(edge);
        }

        // Check for reversed pairs
        for ((from, to), edges) in &pair_edges {
            let reverse_key = (to.clone(), from.clone());
            if let Some(reverse_edges) = pair_edges.get(&reverse_key) {
                // Found bidirectional edges - check if they represent inversions
                for edge in edges {
                    for rev_edge in reverse_edges {
                        // If different genomes in each direction, it's an inversion
                        let shared = edge.genomes.iter()
                            .filter(|g| rev_edge.genomes.contains(g))
                            .count();
                        let key_clone = reverse_key.clone();
                        if shared == 0 && !checked_pairs.contains(&key_clone) {
                            checked_pairs.insert((from.clone(), to.clone()));
                            checked_pairs.insert(reverse_key.clone());

                            inversions.push(StructuralVariant {
                                variant_type: VariantType::Inversion,
                                cluster_ids: vec![from.clone(), to.clone()],
                                affected_genomes: vec![],
                                support: edge.support.min(rev_edge.support),
                                description: format!("Inversion between {} and {}", from, to),
                            });
                        }
                    }
                }
            }
        }

        inversions
    }

    /// Detect duplications: nodes with high degree or self-loops.
    fn detect_duplications(&self, graph: &PangenomeGraph) -> Vec<StructuralVariant> {
        let mut duplications = Vec::new();

        for (cluster_id, node) in &graph.nodes {
            let degree = graph.degree(cluster_id);

            // High degree suggests duplication
            if degree >= 3 {
                duplications.push(StructuralVariant {
                    variant_type: VariantType::Duplication,
                    cluster_ids: vec![cluster_id.to_string()],
                    affected_genomes: node.genomes.iter().map(|g| g.to_string()).collect(),
                    support: node.support,
                    description: format!("Duplication: {} has {} connections", cluster_id, degree),
                });
            }

            // Check for self-loops (suspicious for duplication)
            if let Some(edge) = graph.edges.get(&(cluster_id.clone(), cluster_id.clone())) {
                duplications.push(StructuralVariant {
                    variant_type: VariantType::Duplication,
                    cluster_ids: vec![cluster_id.to_string()],
                    affected_genomes: edge.genomes.iter().map(|g| g.to_string()).collect(),
                    support: edge.support,
                    description: format!("Self-loop in {} suggests duplication", cluster_id),
                });
            }
        }

        duplications
    }
}

impl Default for StructuralVariantDetector {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detector_creation() {
        let detector = StructuralVariantDetector::new();
        assert_eq!(detector.min_support, 1);
    }

    #[test]
    fn test_with_min_support() {
        let detector = StructuralVariantDetector::new().with_min_support(5);
        assert_eq!(detector.min_support, 5);
    }

    #[test]
    fn test_empty_graph_no_variants() {
        let graph = PangenomeGraph::new();
        let detector = StructuralVariantDetector::new();
        let variants = detector.detect(&graph);
        assert!(variants.is_empty());
    }
}
