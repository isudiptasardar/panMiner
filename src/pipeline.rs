//! Main pipeline orchestration for PanMiner.

use rayon::prelude::*;

use crate::config::PanminerConfig;
use crate::clustering::{Clusterer, CpuClusterer, MMseqsRunner};
use crate::correction::{ContaminationRemover, FragmentMerger};
use crate::error::{Error, Result};
use crate::graph::{
    BitPackedMatrix, ConcurrentGraph, Gene, GeneCluster, Node, GenomeId, GraphBuilder, PangenomeGraph,
};
use crate::io::GffParser;
use crate::output::{OutputPaths, OutputWriter};

/// Main pipeline for PanMiner pangenome analysis.
pub struct PanminerPipeline {
    config: PanminerConfig,
}

impl PanminerPipeline {
    /// Create a new pipeline with the given configuration.
    pub fn new(config: PanminerConfig) -> Self {
        Self { config }
    }

    /// Run the full pangenome analysis pipeline.
    pub fn run(&self) -> Result<OutputPaths> {
        self.config.validate()?;

        // Configure rayon thread pool
        let threads = self.config.effective_threads();
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok(); // Ignore error if already initialized

        tracing::info!("Using {} threads", threads);

        // Phase 1: Parse input files
        tracing::info!("Phase 1: Parsing {} input files", self.config.input_files.len());
        let (genes, genome_ids) = self.parse_inputs()?;
        tracing::info!("Parsed {} genes from {} genomes", genes.len(), genome_ids.len());

        if genes.is_empty() {
            return Err(Error::InvalidInput("No genes found in input files".to_string()));
        }

        // Phase 2: Cluster genes
        tracing::info!("Phase 2: Clustering genes");
        let clusters = self.cluster_genes(&genes)?;
        tracing::info!("Found {} clusters", clusters.len());

        // Phase 3: Build graph
        let concurrent_graph = if self.config.chunk_size > 0 && self.config.input_files.len() > self.config.chunk_size {
            tracing::info!("Phase 3: Building pangenome graph (chunked streaming)");
            let streaming = crate::io::StreamingPipeline::new(self.config.clone());
            let chunk_files = streaming.process_chunks_with_clusters(&self.config.input_files, &clusters)?;

            let graph = ConcurrentGraph::with_capacity(clusters.len());
            
            let gene_to_genome: std::collections::HashMap<String, GenomeId> = genes
                .iter()
                .map(|g| (g.id.to_string(), g.genome_id.clone()))
                .collect();

            clusters.par_iter().for_each(|cluster| {
                let mut node = Node::from_cluster(cluster);
                for gene_id in &cluster.genes {
                    if let Some(genome_id) = gene_to_genome.get(gene_id.as_str()) {
                        node.genomes.insert(genome_id.clone());
                    }
                }
                graph.add_node(node);
            });

            for chunk_file in &chunk_files {
                let partial = streaming.load_intermediate(chunk_file)?;
                graph.merge_from(vec![partial]);
            }
            graph
        } else {
            tracing::info!("Phase 3: Building pangenome graph (in-memory)");
            self.build_graph(&clusters, &genes)
        };

        tracing::info!(
            "Graph: {} nodes, {} edges",
            concurrent_graph.node_count(),
            concurrent_graph.edge_count()
        );

        // Phase 4: Error correction
        tracing::info!("Phase 4: Running error correction");
        self.run_corrections(&concurrent_graph, genome_ids.len())?;
        tracing::info!(
            "After correction: {} nodes, {} edges",
            concurrent_graph.node_count(),
            concurrent_graph.edge_count()
        );

        // Phase 5: Build presence/absence matrix
        tracing::info!("Phase 5: Building presence/absence matrix");
        let graph = concurrent_graph.to_standard();
        let matrix = self.build_matrix(&graph, &genome_ids);

        // Phase 6: Generate outputs
        tracing::info!("Phase 6: Generating outputs");
        let writer = OutputWriter::new(&self.config);
        let paths = writer.write_all(&graph, &matrix)?;

        tracing::info!("Pipeline complete");
        Ok(paths)
    }

    /// Parse all input GFF3 files in parallel.
    fn parse_inputs(&self) -> Result<(Vec<Gene>, Vec<GenomeId>)> {
        let results: Vec<(Vec<Gene>, GenomeId)> = self.config.input_files
            .par_iter()
            .map(|path| {
                let genome_id = GenomeId::new(
                    path.file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("unknown")
                );

                let genes = GffParser::open(path, genome_id.clone())
                    .and_then(|p: crate::io::GffParser| p.parse_genes())
                    .unwrap_or_default();

                if genes.is_empty() {
                    tracing::warn!("No genes found in {:?}", path);
                }

                (genes, genome_id)
            })
            .collect();

        let mut all_genes = Vec::new();
        let mut genome_ids = Vec::new();

        for (genes, genome_id) in results {
            all_genes.extend(genes);
            genome_ids.push(genome_id);
        }

        Ok((all_genes, genome_ids))
    }

    /// Cluster genes using MMseqs2-GPU or CPU fallback.
    fn cluster_genes(&self, genes: &[Gene]) -> Result<Vec<GeneCluster>> {
        // Try MMseqs2 first
        if self.config.enable_mmseqs && !self.config.force_cpu {
            if let Some(runner) = MMseqsRunner::detect() {
                tracing::info!("Using {} for clustering", runner.name());
                return runner.cluster(genes, self.config.cluster_identity);
            }
            tracing::info!("MMseqs2 not found, falling back to CPU clustering");
        }

        // CPU fallback
        let clusterer = CpuClusterer::new(self.config.effective_threads());
        tracing::info!("Using {} for clustering", clusterer.name());
        clusterer.cluster(genes, self.config.cluster_identity)
    }

    /// Build the pangenome graph from clusters and genes.
    fn build_graph(&self, clusters: &[GeneCluster], genes: &[Gene]) -> ConcurrentGraph {
        let builder = GraphBuilder::new()
            .with_min_support(self.config.min_support);

        builder.build_concurrent(clusters, genes)
    }

    /// Run error correction on the graph.
    fn run_corrections(&self, graph: &ConcurrentGraph, num_genomes: usize) -> Result<()> {
        // Contamination removal
        let remover = ContaminationRemover::from_mode(&self.config.mode, num_genomes);
        remover.remove(graph)?;

        // Fragment merging (uses empty sequences for now)
        // TODO: pass actual cluster sequences when available
        let merger = FragmentMerger::new()
            .with_collapse_threshold(self.config.collapse_threshold);

        let empty_sequences = std::collections::HashMap::new();
        merger.correct_mistranslations(graph, &empty_sequences)?;
        merger.collapse_gene_families(graph, &empty_sequences)?;

        Ok(())
    }

    /// Build presence/absence matrix from the graph.
    fn build_matrix(&self, graph: &PangenomeGraph, genome_ids: &[GenomeId]) -> BitPackedMatrix {
        let num_genomes = genome_ids.len();
        let cluster_ids: Vec<_> = graph.nodes.keys().collect();
        let num_clusters = cluster_ids.len();

        let mut matrix = BitPackedMatrix::new(num_genomes, num_clusters);

        // Set genome names
        matrix.set_genome_names(genome_ids.iter().map(|id| id.to_string()).collect());
        matrix.set_cluster_ids(cluster_ids.iter().map(|id| id.to_string()).collect());

        // Create a fast lookup for genome index
        let genome_to_idx: std::collections::HashMap<&GenomeId, usize> = genome_ids
            .iter()
            .enumerate()
            .map(|(idx, id)| (id, idx))
            .collect();

        // Populate matrix accurately using node.genomes
        for (cluster_idx, cluster_id) in cluster_ids.iter().enumerate() {
            if let Some(node) = graph.nodes.get(*cluster_id) {
                for genome_id in &node.genomes {
                    if let Some(&genome_idx) = genome_to_idx.get(genome_id) {
                        matrix.set(genome_idx, cluster_idx, true);
                    }
                }
            }
        }

        matrix
    }
}

#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn test_pipeline_creation() {
        let config = PanminerConfig::default();
        let _pipeline = PanminerPipeline::new(config);
        // Just verify it doesn't panic
    }

    #[test]
    fn test_pipeline_validation() {
        let config = PanminerConfig::default(); // No input files
        let pipeline = PanminerPipeline::new(config);
        assert!(pipeline.run().is_err());
    }
}

#[cfg(test)]
mod streaming_tests {
    use super::*;


    #[test]
    fn test_pipeline_chunked_streaming() {
        let temp_dir = tempfile::tempdir().unwrap();
        let dir = temp_dir.path();
        let gff1 = dir.join("seq1.gff");
        std::fs::write(&gff1, "##gff-version 3\nseq1\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene1;product=test\n##FASTA\n>seq1\nATCGATCGATCGATCG\n").unwrap();
        let gff2 = dir.join("seq2.gff");
        std::fs::write(&gff2, "##gff-version 3\nseq2\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene2;product=test\n##FASTA\n>seq2\nATCGATCGATCGATCG\n").unwrap();
        let gff3 = dir.join("seq3.gff");
        std::fs::write(&gff3, "##gff-version 3\nseq3\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene3;product=test\n##FASTA\n>seq3\nATCGATCGATCGATCG\n").unwrap();

        let config = PanminerConfig::default()
            .with_input_files(vec![
                gff1.clone(),
                gff2.clone(),
                gff3.clone(),
            ])
            .with_output_dir(temp_dir.path().to_path_buf())
            .with_temp_dir(temp_dir.path().to_path_buf())
            .with_chunk_size(2); // 3 files, chunk size 2 -> 2 chunks

        let pipeline = PanminerPipeline::new(config);
        let paths = pipeline.run().expect("Pipeline should run successfully with chunks");

        // The matrix should have 3 genomes
        let matrix_content = std::fs::read_to_string(paths.matrix.as_ref().unwrap()).unwrap();
        assert!(matrix_content.contains("seq1"));
        assert!(matrix_content.contains("seq2"));
        assert!(matrix_content.contains("seq3"));
        
        // Verify intermediate files were created
        let intermediate_dir = temp_dir.path().join("panminer_intermediate");
        assert!(intermediate_dir.exists(), "Intermediate directory should exist");
        let mut chunk_files = 0;
        for entry in std::fs::read_dir(intermediate_dir).unwrap() {
            let entry = entry.unwrap();
            if entry.path().extension().and_then(|s| s.to_str()) == Some("zst") {
                chunk_files += 1;
            }
        }
        assert_eq!(chunk_files, 2, "Should have created 2 chunk files");
    }
}
