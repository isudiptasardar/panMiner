//! Graph construction from gene clusters.

use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

use super::concurrent::ConcurrentGraph;
use super::types::{ClusterId, Gene, GeneCluster, GeneId, GenomeId, Node, PangenomeGraph};


/// Builder for constructing pangenome graphs from gene clusters.
///
/// This builder creates the graph structure from clustered genes,
/// handling parallel construction and edge creation.
pub struct GraphBuilder {
    /// Minimum support threshold for keeping nodes
    min_support: usize,
    /// Full contig DNA from GFF FASTA section (keyed by genome, contig name)
    contig_dna: Option<HashMap<(GenomeId, String), Vec<u8>>>,
}

impl GraphBuilder {
    /// Create a new graph builder with default settings.
    pub fn new() -> Self {
        Self { min_support: 1, contig_dna: None }
    }

    /// Set minimum support threshold.
    pub fn with_min_support(mut self, min_support: usize) -> Self {
        self.min_support = min_support;
        self
    }

    /// Set full contig DNA from GFF FASTA section.
    ///
    /// When provided, nodes will store full contig DNA (including intergenic
    /// regions) instead of concatenated gene sequences. This is needed for
    /// missing gene recovery which searches flanking regions.
    pub fn with_contig_dna(mut self, contig_dna: &HashMap<(GenomeId, String), Vec<u8>>) -> Self {
        self.contig_dna = Some(contig_dna.clone());
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

        // Build gene data map for populating gene_members on nodes
        let gene_data_map: HashMap<GeneId, Gene> = genes.iter()
            .map(|g| (g.id.clone(), g.clone()))
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

        // Identify contig-end genes (first and last gene on each contig)
        let mut contig_end_gene_ids: HashSet<String> = HashSet::new();
        for ((_genome, _contig), contig_genes) in &genes_by_contig {
            if contig_genes.len() == 1 {
                // Single gene on contig → both start and end
                contig_end_gene_ids.insert(contig_genes[0].id.to_string());
            } else {
                // Multiple genes: sort by position, mark first and last
                let mut sorted: Vec<_> = contig_genes.iter().collect();
                sorted.sort_by_key(|g| g.start);
                contig_end_gene_ids.insert(sorted.first().unwrap().id.to_string());
                contig_end_gene_ids.insert(sorted.last().unwrap().id.to_string());
            }
        }

        // Extract contig sequences by concatenating gene sequences
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

        // Build a mapping: cluster_id → set of (genome, contig) pairs it appears on
        let cluster_to_contigs: HashMap<ClusterId, HashSet<(GenomeId, String)>> = clusters
            .iter()
            .map(|cluster| {
                let contigs: HashSet<(GenomeId, String)> = cluster.genes.iter()
                    .filter_map(|gid| {
                        let genome = gene_to_genome.get(gid.as_str())?;
                        let gene = gene_data_map.get(gid)?;
                        Some((genome.clone(), gene.contig.clone()))
                    })
                    .collect();
                (cluster.id.clone(), contigs)
            })
            .collect();

        // Add nodes from clusters (parallel)
        clusters.par_iter().for_each(|cluster| {
            let mut node = Node::from_cluster_with_genes(cluster, &gene_data_map);
            for gene_id in &cluster.genes {
                if let Some(genome_id) = gene_to_genome.get(gene_id.as_str()) {
                    node.genomes.insert(genome_id.clone());
                }
            }
            // Track which genomes have this cluster at a contig boundary
            for gene_id in &cluster.genes {
                if contig_end_gene_ids.contains(gene_id.as_str()) {
                    if let Some(genome_id) = gene_to_genome.get(gene_id.as_str()) {
                        node.contig_end_genomes.insert(genome_id.clone());
                    }
                }
            }
            // Populate contig_sequences: prefer full contig DNA from GFF FASTA,
            // fall back to concatenated gene sequences
            if let Some(contigs) = cluster_to_contigs.get(&cluster.id) {
                for (genome, contig_name) in contigs {
                    // Try full contig DNA first (from GFF FASTA section)
                    if let Some(full_dna) = self.contig_dna.as_ref()
                        .and_then(|cd| cd.get(&(genome.clone(), contig_name.clone())))
                    {
                        node.add_contig_sequence(contig_name.clone(), full_dna.clone());
                    } else if let Some(seq) = contig_sequences.get(&(genome.clone(), contig_name.clone())) {
                        // Fall back to concatenated gene sequences
                        node.add_contig_sequence(contig_name.clone(), seq.clone());
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
        let mut graph = concurrent.to_standard();

        // Populate gene lookup for output writers
        for gene in genes {
            graph.gene_lookup.insert(gene.id.clone(), gene.clone());
        }

        graph
    }

    /// Build a graph from pre-extracted adjacencies.
    ///
    /// This is used when loading from intermediate files.
    /// The adjacencies tuple is (genome_id, from_cluster, to_cluster).
    pub fn build_from_adjacencies(
        &self,
        adjacencies: &[(String, String, String)], // (genome_id, from_cluster, to_cluster)
        _genomes: &[GenomeId],
    ) -> ConcurrentGraph {
        let graph = ConcurrentGraph::new();

        // Add edges from adjacencies
        adjacencies.par_iter().for_each(|(genome_id, from, to)| {
            let from_cluster = ClusterId::new(from);
            let to_cluster = ClusterId::new(to);
            let genome = GenomeId::new(genome_id);

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

    #[test]
    fn test_contig_end_marking() {
        use super::super::types::{GeneId, Strand};

        let mut cluster1 = GeneCluster::new("c1");
        cluster1.support = 5;
        cluster1.add_gene(GeneId::new("g1"));

        let mut cluster2 = GeneCluster::new("c2");
        cluster2.support = 3;
        cluster2.add_gene(GeneId::new("g2"));

        let mut cluster3 = GeneCluster::new("c3");
        cluster3.support = 2;
        cluster3.add_gene(GeneId::new("g3"));

        let clusters = vec![cluster1, cluster2, cluster3];

        // Two contigs: contig1 has g1+g2, contig2 has g3 alone
        let gene1 = Gene {
            id: GeneId::new("g1"), sequence: b"ATCG".to_vec(),
            genome_id: GenomeId::new("genome1"), contig: "contig1".to_string(),
            start: 1, end: 4, strand: Strand::Plus, annotation: None,
        };
        let gene2 = Gene {
            id: GeneId::new("g2"), sequence: b"GCTA".to_vec(),
            genome_id: GenomeId::new("genome1"), contig: "contig1".to_string(),
            start: 10, end: 13, strand: Strand::Plus, annotation: None,
        };
        let gene3 = Gene {
            id: GeneId::new("g3"), sequence: b"TTTT".to_vec(),
            genome_id: GenomeId::new("genome1"), contig: "contig2".to_string(),
            start: 1, end: 4, strand: Strand::Plus, annotation: None,
        };

        let genes = vec![gene1, gene2, gene3];
        let graph = GraphBuilder::new().build_concurrent(&clusters, &genes);

        // g1 and g2 are on contig1 — g1 is first (start=1), g2 is last (start=10)
        // Both should be marked as contig-end
        let node_c1 = graph.nodes.get(&ClusterId::new("c1")).unwrap();
        let node_c2 = graph.nodes.get(&ClusterId::new("c2")).unwrap();
        let node_c3 = graph.nodes.get(&ClusterId::new("c3")).unwrap();

        assert!(node_c1.contig_end_genomes.contains(&GenomeId::new("genome1")), "first gene on contig should be marked as contig end");
        assert!(node_c2.contig_end_genomes.contains(&GenomeId::new("genome1")), "last gene on contig should be marked as contig end");
        assert!(node_c3.contig_end_genomes.contains(&GenomeId::new("genome1")), "sole gene on contig should be marked as contig end");
    }

    #[test]
    fn test_contig_sequences_populated() {
        use super::super::types::{GeneId, Strand};

        let mut cluster1 = GeneCluster::new("c1");
        cluster1.support = 1;
        cluster1.add_gene(GeneId::new("g1"));

        let clusters = vec![cluster1];

        let gene1 = Gene {
            id: GeneId::new("g1"), sequence: b"ATCG".to_vec(),
            genome_id: GenomeId::new("genome1"), contig: "contig1".to_string(),
            start: 1, end: 4, strand: Strand::Plus, annotation: None,
        };

        let genes = vec![gene1];
        let graph = GraphBuilder::new().build_concurrent(&clusters, &genes);

        let node_c1 = graph.nodes.get(&ClusterId::new("c1")).unwrap();
        assert!(!node_c1.contig_sequences.is_empty(),
            "contig_sequences should be populated from gene sequences");
    }
}