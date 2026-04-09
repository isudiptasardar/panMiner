//! Graph construction from gene clusters.

use rayon::prelude::*;
use std::collections::HashMap;

use super::concurrent::ConcurrentGraph;
use super::types::{ClusterId, Gene, GeneCluster, GenomeId, Node, PangenomeGraph};


/// Builder for constructing pangenome graphs from gene clusters.
///
/// This builder creates the graph structure from clustered genes,
/// handling parallel construction and edge creation.
pub struct GraphBuilder {
    /// Minimum support threshold for keeping nodes
    min_support: usize,
}

impl GraphBuilder {
    /// Create a new graph builder with default settings.
    pub fn new() -> Self {
        Self { min_support: 1 }
    }

    /// Set minimum support threshold.
    pub fn with_min_support(mut self, min_support: usize) -> Self {
        self.min_support = min_support;
        self
    }

    /// Build a concurrent graph from gene clusters and genes.
    ///
    /// This method constructs the graph in parallel using DashMap,
    /// which allows lock-free concurrent updates.
    pub fn build_concurrent(
        &self,
        clusters: &[GeneCluster],
        genes: &[Gene],
    ) -> ConcurrentGraph {
        let graph = ConcurrentGraph::with_capacity(clusters.len());

        // Create gene to cluster mapping and gene to genome mapping
        let gene_to_cluster: HashMap<String, ClusterId> = clusters
            .iter()
            .flat_map(|cluster| {
                cluster.genes.iter().map(|gene_id| {
                    (gene_id.to_string(), cluster.id.clone())
                })
            })
            .collect();

        let gene_to_genome: HashMap<String, GenomeId> = genes
            .iter()
            .map(|g| (g.id.to_string(), g.genome_id.clone()))
            .collect();

        // Group genes by genome and contig for contig sequence collection
        let genes_by_contig: HashMap<(GenomeId, String), Vec<&Gene>> = genes
            .iter()
            .fold(HashMap::new(), |mut acc, gene| {
                let key = (gene.genome_id.clone(), gene.contig.clone());
                acc.entry(key).or_default().push(gene);
                acc
            });

        // Count genes per (genome, contig) for contig-end marking
        let contig_gene_count: HashMap<(GenomeId, String), usize> = genes_by_contig
            .iter()
            .map(|((genome, contig), genes)| ((genome.clone(), contig.clone()), genes.len()))
            .collect();

        // Extract contig sequences
        let contig_sequences: HashMap<(GenomeId, String), Vec<u8>> = genes_by_contig
            .iter()
            .map(|((genome, contig), genes)| {
                let mut seq = Vec::new();
                for g in genes {
                    seq.extend(&g.sequence);
                }
                ((genome.clone(), contig.clone()), seq)
            })
            .collect();

        // Add nodes from clusters (parallel)
        clusters.par_iter().for_each(|cluster| {
            let mut node = Node::from_cluster(cluster);
            for gene_id in &cluster.genes {
                if let Some(genome_id) = gene_to_genome.get(gene_id.as_str()) {
                    node.genomes.insert(genome_id.clone());
                }
            }
            // Add contig sequences if available and mark contig ends
            for ((genome, contig), seq) in &contig_sequences {
                if node.genomes.contains(genome) {
                    node.add_contig_sequence(contig.clone(), seq.clone());
                    // Mark as contig end if this node represents the only gene on its contig
                    if let Some(genome_id) = node.genomes.iter().next() {
                        if contig_gene_count.get(&(genome_id.clone(), contig.clone())).map(|c| *c == 1).unwrap_or(false) {
                            node.is_contig_end = true;
                        }
                    }
                }
            }
            graph.add_node(node);
        });

        // Build edges from adjacencies (parallel)
        genes_by_contig
            .par_iter()
            .for_each(|((_genome, _contig), contig_genes)| {
                // Sort by position
                let mut sorted: Vec<_> = contig_genes.iter().collect();
                sorted.sort_by_key(|g| g.start);

                // Create edges between adjacent genes
                for window in sorted.windows(2) {
                    let gene1 = window[0];
                    let gene2 = window[1];

                    // Find clusters for each gene
                    if let (Some(c1), Some(c2)) = (
                        gene_to_cluster.get(&gene1.id.to_string()),
                        gene_to_cluster.get(&gene2.id.to_string()),
                    ) {
                        graph.add_edge_genome(c1.clone(), c2.clone(), gene1.genome_id.clone());
                    }
                }
            });

        graph
    }

    /// Build a standard graph from gene clusters and genes.
    ///
    /// This is a convenience method that builds a concurrent graph
    /// and converts it to the standard representation.
    pub fn build(
        &self,
        clusters: &[GeneCluster],
        genes: &[Gene],
    ) -> PangenomeGraph {
        let concurrent = self.build_concurrent(clusters, genes);
        concurrent.to_standard()
    }

    /// Build a graph from pre-extracted adjacencies.
    ///
    /// This is used when loading from intermediate files.
    pub fn build_from_adjacencies(
        &self,
        adjacencies: &[(String, String, String)], // (contig, from_cluster, to_cluster)
        _genomes: &[GenomeId],
    ) -> ConcurrentGraph {
        let graph = ConcurrentGraph::new();

        // Add edges from adjacencies
        adjacencies.par_iter().for_each(|(contig, from, to)| {
            let from_cluster = ClusterId::new(from);
            let to_cluster = ClusterId::new(to);
            let genome = GenomeId::new(contig); // Contig encodes genome

            graph.add_edge_genome(from_cluster, to_cluster, genome);
        });

        graph
    }

    /// Filter low-support nodes from the graph.
    ///
    /// This removes nodes that appear in fewer than `min_support` genomes
    /// and have degree <= 1 (only one connected edge).
    pub fn filter_low_support(&self, graph: &ConcurrentGraph, threshold: usize) -> Vec<ClusterId> {
        graph.find_low_support_nodes(threshold)
    }
}

impl Default for GraphBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_builder_construction() {
        let builder = GraphBuilder::new();

        // Create test clusters
        let mut cluster1 = GeneCluster::new("c1");
        cluster1.support = 5;
        cluster1.add_gene(GeneId::new("g1"));

        let mut cluster2 = GeneCluster::new("c2");
        cluster2.support = 3;
        cluster2.add_gene(GeneId::new("g2"));

        let clusters = vec![cluster1, cluster2];

        // Create test genes
        use super::super::types::{GeneId, Strand};

        let gene1 = Gene {
            id: GeneId::new("g1"),
            sequence: b"ATCG".to_vec(),
            genome_id: GenomeId::new("genome1"),
            contig: "contig1".to_string(),
            start: 1,
            end: 4,
            strand: Strand::Plus,
            annotation: None,
        };

        let gene2 = Gene {
            id: GeneId::new("g2"),
            sequence: b"GCTA".to_vec(),
            genome_id: GenomeId::new("genome1"),
            contig: "contig1".to_string(),
            start: 10,
            end: 13,
            strand: Strand::Plus,
            annotation: None,
        };

        let genes = vec![gene1, gene2];

        let graph = builder.build(&clusters, &genes);

        assert_eq!(graph.node_count(), 2);
    }
}