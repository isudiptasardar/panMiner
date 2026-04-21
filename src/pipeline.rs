//! Main pipeline orchestration for PanMiner.

use rayon::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;

use crate::config::{CorrectionMode, PanminerConfig, PipelineMode};
use crate::clustering::{Clusterer, CpuClusterer, MMseqsRunner};
use crate::correction::{ContaminationRemover, ContigEndPruner, FragmentMerger, MissingGeneRecoverer, MisassemblyEdgeCleaner, ParalogResolver, DistanceCache};
use crate::error::{Error, Result};
use crate::graph::{
    BitPackedMatrix, ConcurrentGraph, Gene, GeneCluster, Node, GenomeId, GraphBuilder, PangenomeGraph,
};
#[cfg(feature = "prodigal")]
use crate::graph::GeneId;
use crate::io::{GffParser, CheckmQcRunner, QcRunner, GenomeQC, BaktaRunner, is_genbank_file};
use crate::output::{OutputPaths, OutputWriter};
use crate::output::qc_stats::{write_qc_stats, write_qc_summary};
use crate::gwas::{PyseerRunner, GWASRunner};
use std::fs;

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

        // Dispatch based on pipeline mode
        if self.config.pipeline_mode == PipelineMode::Dbg {
            return self.run_dbg_mode();
        }
        if self.config.pipeline_mode == PipelineMode::Prodigal {
            return self.run_prodigal_mode();
        }

        // Phase 0: Pre-processing QC (optional)
        let qc_results = if self.config.enable_qc {
            tracing::info!("Phase 0: Running pre-processing QC");
            self.run_qc()?
        } else {
            tracing::info!("Phase 0: Skipped (QC disabled)");
            vec![]
        };

        // Filter out genomes that failed QC
        let input_files = if !qc_results.is_empty() {
            if qc_results.len() != self.config.input_files.len() {
                tracing::warn!(
                    "QC results count ({}) doesn't match input files count ({}). Some QC runs may have failed.",
                    qc_results.len(),
                    self.config.input_files.len()
                );
            }
            let passed_files: Vec<PathBuf> = self.config.input_files.iter()
                .enumerate()
                .filter(|(i, _)| {
                    qc_results.get(*i).map(|qc| qc.passed).unwrap_or(true)
                })
                .map(|(_, path)| path.clone())
                .collect();

            let removed = self.config.input_files.len() - passed_files.len();
            if removed > 0 {
                tracing::info!(
                    "QC filtering: {} genomes passed, {} removed",
                    passed_files.len(),
                    removed
                );
            } else {
                tracing::info!("QC filtering: all {} genomes passed", passed_files.len());
            }

            if passed_files.is_empty() {
                return Err(Error::InvalidInput(
                    "No genomes passed QC filtering. Use --no-qc to disable.".to_string()
                ));
            }

            passed_files
        } else {
            self.config.input_files.clone()
        };
        let input_files = if self.config.reannotate {
            tracing::info!("Phase 0.5: Re-annotating inputs with Bakta");
            self.reannotate_inputs(&input_files)?
        } else {
            input_files
        };

        // Phase 1: Parse input files
        tracing::info!("Phase 1: Parsing {} input files", input_files.len());
        let (genes, genome_ids, contig_dna) = self.parse_inputs_from(&input_files)?;
        tracing::info!("Parsed {} genes from {} genomes", genes.len(), genome_ids.len());

        if genes.is_empty() {
            return Err(Error::InvalidInput("No genes found in input files".to_string()));
        }

        // Phase 2: Cluster genes
        tracing::info!("Phase 2: Clustering genes");
        let mut clusters = self.cluster_genes(&genes)?;
        tracing::info!("Found {} clusters", clusters.len());

        // Mark paralog clusters (same genome appears multiple times in a cluster)
        ParalogResolver::mark_paralog_clusters(&mut clusters, &genes);
        let paralog_count = clusters.iter().filter(|c| c.is_paralog).count();
        if paralog_count > 0 {
            tracing::info!("Marked {} clusters as containing paralogs", paralog_count);
        }

        // Phase 3: Build graph
        let concurrent_graph = if self.config.chunk_size > 0 && input_files.len() > self.config.chunk_size {
            tracing::info!("Phase 3: Building pangenome graph (chunked streaming)");
            let streaming = crate::io::StreamingPipeline::new(self.config.clone());
            let chunk_files = streaming.process_chunks_with_clusters(&input_files, &clusters)?;

            let graph = ConcurrentGraph::with_capacity(clusters.len());

            let gene_to_genome: std::collections::HashMap<String, GenomeId> = genes
                .iter()
                .map(|g| (g.id.to_string(), g.genome_id.clone()))
                .collect();

            // Build gene data map for O(1) lookup (same as GraphBuilder)
            let gene_data_map: std::collections::HashMap<crate::graph::GeneId, Gene> = genes
                .iter()
                .map(|g| (g.id.clone(), g.clone()))
                .collect();

            clusters.par_iter().for_each(|cluster| {
                let mut node = Node::from_cluster(cluster);
                for gene_id in &cluster.genes {
                    if let Some(genome_id) = gene_to_genome.get(gene_id.as_str()) {
                        node.genomes.insert(genome_id.clone());
                    }
                }
                // Add contig DNA to nodes using O(1) gene data lookup
                for gene_id in &cluster.genes {
                    if let Some(genome_id) = gene_to_genome.get(gene_id.as_str()) {
                        if let Some(gene) = gene_data_map.get(gene_id) {
                            let key = (genome_id.clone(), gene.contig.clone());
                            if let Some(full_dna) = contig_dna.get(&key) {
                                node.add_contig_sequence(gene.contig.clone(), full_dna.clone());
                            }
                        }
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
            self.build_graph(&clusters, &genes, &contig_dna)
        };

        tracing::info!(
            "Graph: {} nodes, {} edges",
            concurrent_graph.node_count(),
            concurrent_graph.edge_count()
        );

        // Write pre-filtered graph for debugging/comparison
        fs::create_dir_all(&self.config.output_dir)?;
        {
            let pre_graph = concurrent_graph.to_standard();
            let pre_graph_path = self.config.output_dir.join("pre_filt_graph.gml");
            if let Err(e) = crate::output::GmlWriter::write(&pre_graph, &pre_graph_path) {
                tracing::warn!("Failed to write pre-filtered graph: {}", e);
            } else {
                tracing::info!("Wrote pre_filt_graph.gml");
            }
        }

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
        let mut graph = concurrent_graph.to_standard();
        let matrix = self.build_matrix(&graph, &genome_ids);

        // Phase 5.5: Detect highly variable genes
        tracing::info!("Phase 5.5: Detecting highly variable genes");
        let hv_detector = crate::graph::HighlyVariableDetector::new();
        let hv_result = hv_detector.detect(&graph);
        for cluster_id in &hv_result.highly_variable {
            if let Some(node) = graph.nodes.get_mut(cluster_id) {
                node.is_highly_variable = true;
            }
        }
        tracing::info!(
            "Highly variable detection: {} cycles found, {} merged sets, {} genes flagged",
            hv_result.cycles_found, hv_result.merged_sets, hv_result.highly_variable.len()
        );

        // Phase 6: Generate outputs
        tracing::info!("Phase 6: Generating outputs");
        let writer = OutputWriter::new(&self.config);

        // Build gene_members map for Roary CSV output
        let gene_members = graph.build_gene_members_map();

        let mut paths = writer.write_all(&graph, &matrix, &gene_members)?;

        // Track pre-filtered graph path (written before corrections)
        let pre_filt_path = self.config.output_dir.join("pre_filt_graph.gml");
        if pre_filt_path.exists() {
            paths.pre_filt_graph = Some(pre_filt_path);
        }

        // Phase 7: GWAS analysis (optional)
        if self.config.run_gwas {
            tracing::info!("Phase 7: Running GWAS analysis");
            self.run_gwas(&graph, &matrix)?;
        }

        // Clean up Bakta temporary files
        if self.config.reannotate && !self.config.keep_bakta_output {
            let bakta_tmp = self.config.output_dir.join("bakta_tmp");
            if bakta_tmp.exists() {
                if let Err(e) = std::fs::remove_dir_all(&bakta_tmp) {
                    tracing::warn!("Failed to clean up Bakta temp files: {}", e);
                } else {
                    tracing::info!("Cleaned up Bakta temporary files");
                }
            }
        }

        tracing::info!("Pipeline complete");
        Ok(paths)
    }

    /// Run the cDBG-based pipeline (GGCAT + ggCaller).
    ///
    /// Flow:
    /// 1. Build colored cDBG with GGCAT (feature-gated)
    /// 2. Call genes with ggCaller (subprocess)
    /// 3. Parse ggCaller GFF output into PanMiner Gene structs
    /// 4. Build pangenome graph (reuses existing code)
    /// 5. Run corrections (reuses existing code)
    /// 6. Generate outputs (reuses existing code)
    #[allow(unused_variables)]
    fn run_dbg_mode(&self) -> Result<OutputPaths> {
        tracing::info!("Running cDBG-based pipeline (mode=dbg)");

        // --- Phase 0: Pre-processing QC (optional, reuses existing) ---
        let qc_results = if self.config.enable_qc {
            tracing::info!("Phase 0: Running pre-processing QC");
            self.run_qc()?
        } else {
            tracing::info!("Phase 0: Skipped (QC disabled)");
            vec![]
        };

        // --- Phase 1: Build colored cDBG with GGCAT ---
        let genomes: Vec<(PathBuf, String)> = self.config.input_files
            .iter()
            .filter_map(|p| {
                let name = p.file_stem()
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_else(|| "unknown".to_string());
                Some((p.clone(), name))
            })
            .collect();

        if genomes.is_empty() {
            return Err(Error::NoGenomes);
        }

        #[cfg(feature = "dbg")]
        let _cdbg_graph = {
            let cdbg_output_dir = self.config.output_dir.join("cdbg");
            let builder = crate::io::GGCATBuilder::new()
                .with_threads(self.config.effective_threads())
                .with_kmer_size(self.config.kmer_size);
            builder.build_colored_cdbg(&genomes, self.config.kmer_size, &cdbg_output_dir)?
        };

        #[cfg(not(feature = "dbg"))]
        {
            return Err(Error::FeatureNotEnabled("dbg".to_string()));
        }

        // The following code is only reachable when the dbg feature is enabled.
        // The #[cfg(not(feature = "dbg"))] block above returns early, making
        // this code unreachable without the feature. #[allow] suppresses the
        // warning for the default (no dbg) build.
        #[allow(unreachable_code)]

        // --- Phase 2: Call genes with ggCaller ---
        let ggcaller_runner = crate::io::GGCallerRunner::detect()
            .ok_or_else(Error::ggcaller_not_found)?;

        let ggcaller_output_dir = self.config.output_dir.join("ggcaller_output");
        let ggcaller_output = ggcaller_runner.call_genes(
            &self.config.input_files,
            &ggcaller_output_dir,
            self.config.effective_threads(),
        )?;

        // --- Phase 3: Parse ggCaller GFF output ---
        let gff_files = crate::io::GGCallerRunner::parse_gff_paths(&ggcaller_output)?;

        tracing::info!(
            "Parsed {} GFF files from ggCaller output",
            gff_files.len(),
        );

        if gff_files.is_empty() {
            return Err(Error::InvalidInput(
                "ggCaller produced no GFF output files".to_string(),
            ));
        }

        // --- Phase 4: Parse GFF files (reuses existing) ---
        tracing::info!("Phase 4: Parsing {} GFF files from ggCaller", gff_files.len());
        let (genes, genome_ids, contig_dna) = self.parse_inputs_from(&gff_files)?;
        tracing::info!("Parsed {} genes from {} genomes", genes.len(), genome_ids.len());

        if genes.is_empty() {
            return Err(Error::InvalidInput("No genes found in ggCaller GFF output".to_string()));
        }

        // --- Phase 5: Cluster genes (reuses existing) ---
        tracing::info!("Phase 5: Clustering genes");
        let mut clusters = self.cluster_genes(&genes)?;
        tracing::info!("Found {} clusters", clusters.len());

        // Mark paralog clusters
        ParalogResolver::mark_paralog_clusters(&mut clusters, &genes);
        let paralog_count = clusters.iter().filter(|c| c.is_paralog).count();
        if paralog_count > 0 {
            tracing::info!("Marked {} clusters as containing paralogs", paralog_count);
        }

        // --- Phase 6: Build pangenome graph (reuses existing) ---
        tracing::info!("Phase 6: Building pangenome graph");
        let concurrent_graph = self.build_graph(&clusters, &genes, &contig_dna);
        tracing::info!(
            "Graph: {} nodes, {} edges",
            concurrent_graph.node_count(),
            concurrent_graph.edge_count()
        );

        // --- Phase 7: Run corrections (reuses existing) ---
        tracing::info!("Phase 7: Running error correction");
        self.run_corrections(&concurrent_graph, genome_ids.len())?;
        tracing::info!(
            "After correction: {} nodes, {} edges",
            concurrent_graph.node_count(),
            concurrent_graph.edge_count()
        );

        // --- Phase 8: Build presence/absence matrix (reuses existing) ---
        tracing::info!("Phase 8: Building presence/absence matrix");
        let mut graph = concurrent_graph.to_standard();
        let matrix = self.build_matrix(&graph, &genome_ids);

        // --- Phase 8.5: Detect highly variable genes ---
        tracing::info!("Phase 8.5: Detecting highly variable genes");
        let hv_detector = crate::graph::HighlyVariableDetector::new();
        let hv_result = hv_detector.detect(&graph);
        for cluster_id in &hv_result.highly_variable {
            if let Some(node) = graph.nodes.get_mut(cluster_id) {
                node.is_highly_variable = true;
            }
        }
        tracing::info!(
            "Highly variable detection: {} cycles found, {} merged sets, {} genes flagged",
            hv_result.cycles_found, hv_result.merged_sets, hv_result.highly_variable.len()
        );

        // --- Phase 9: Generate outputs (reuses existing) ---
        tracing::info!("Phase 9: Generating outputs");
        let writer = OutputWriter::new(&self.config);

        // Build gene_members map for Roary CSV output
        let gene_members = graph.build_gene_members_map();

        let mut output_paths = writer.write_all(&graph, &matrix, &gene_members)?;

        // Track pre-filtered graph path (if written)
        let pre_filt_path = self.config.output_dir.join("pre_filt_graph.gml");
        if pre_filt_path.exists() {
            output_paths.pre_filt_graph = Some(pre_filt_path);
        }

        // Write QC results if any
        if !qc_results.is_empty() {
            let stats_path = self.config.output_dir.join("qc_stats.csv");
            if let Err(e) = write_qc_stats(&qc_results, &stats_path) {
                tracing::warn!("Failed to write QC stats: {}", e);
            }
        }

        tracing::info!("cDBG pipeline complete. Output: {:?}", output_paths.output_dir);
        Ok(output_paths)
    }

    /// Run Prodigal-based pipeline: call genes on raw FASTA assemblies.
    fn run_prodigal_mode(&self) -> Result<OutputPaths> {
        tracing::info!("Running Prodigal-based pipeline (mode=prodigal)");

        #[cfg(feature = "prodigal")]
        {
            let runner = crate::io::OrphosRunner::detect()
                .ok_or_else(|| Error::ExternalTool(
                    "prodigal not found: install with conda install -c bioconda prodigal".into()
                ))?;

            // Phase 0: Pre-processing QC (optional)
            let qc_results = if self.config.enable_qc {
                tracing::info!("Phase 0: Running pre-processing QC");
                self.run_qc()?
            } else {
                vec![]
            };

            let input_files: Vec<PathBuf> = if !qc_results.is_empty() {
                self.config.input_files.iter()
                    .zip(qc_results.iter())
                    .filter(|(_, qc)| qc.passed)
                    .map(|(p, _)| p.clone())
                    .collect()
            } else {
                self.config.input_files.clone()
            };

            if input_files.is_empty() {
                return Err(Error::NoGenomes);
            }

            // Phase 1: Run Prodigal on each FASTA assembly
            tracing::info!("Phase 1: Running Prodigal gene calling on {} assemblies", input_files.len());
            let mut all_genes = Vec::new();
            let mut genome_ids = Vec::new();
            let mut contig_dna: HashMap<(GenomeId, String), Vec<u8>> = HashMap::new();

            for fasta_path in &input_files {
                let genome_id = fasta_path.file_stem()
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_else(|| "unknown".to_string());
                let gid = GenomeId::new(&genome_id);
                genome_ids.push(gid.clone());

                match runner.predict_genes(fasta_path) {
                    Ok(predicted) => {
                        tracing::info!("Prodigal predicted {} genes for {}", predicted.len(), genome_id);
                        for pg in predicted {
                            let strand = match pg.strand {
                                crate::io::orphos::Strand::Forward => crate::graph::Strand::Plus,
                                crate::io::orphos::Strand::Reverse => crate::graph::Strand::Minus,
                            };
                            let gene = Gene {
                                id: GeneId::new(&pg.gene_id),
                                sequence: pg.sequence,
                                genome_id: gid.clone(),
                                contig: pg.contig,
                                start: pg.start,
                                end: pg.end,
                                strand,
                                annotation: None,
                            };
                            all_genes.push(gene);
                        }
                    }
                    Err(e) => {
                        tracing::warn!("Prodigal failed for {}: {}", genome_id, e);
                    }
                }

                // Read FASTA contig sequences for missing gene recovery
                if let Ok(parser) = crate::io::FastaParser::open(fasta_path) {
                    if let Ok(records) = parser.parse_all() {
                        for (contig_name, seq) in records {
                            contig_dna.insert((gid.clone(), contig_name), seq);
                        }
                    }
                }
            }

            if all_genes.is_empty() {
                return Err(Error::NoGenes(format!(
                    "Prodigal predicted no genes from {} assemblies", input_files.len()
                )));
            }

            tracing::info!("Prodigal produced {} genes from {} genomes", all_genes.len(), genome_ids.len());

            // Phases 2-9: same as GFF pipeline
            let mut clusters = self.cluster_genes(&all_genes)?;
            ParalogResolver::mark_paralog_clusters(&mut clusters, &all_genes);
            let concurrent_graph = self.build_graph(&clusters, &all_genes, &contig_dna);
            self.run_corrections(&concurrent_graph, genome_ids.len())?;

            let mut graph = concurrent_graph.to_standard();
            let matrix = self.build_matrix(&graph, &genome_ids);

            let hv_detector = crate::graph::HighlyVariableDetector::new();
            let hv_result = hv_detector.detect(&graph);
            for cluster_id in &hv_result.highly_variable {
                if let Some(node) = graph.nodes.get_mut(cluster_id) {
                    node.is_highly_variable = true;
                }
            }

            let writer = OutputWriter::new(&self.config);
            let gene_members = graph.build_gene_members_map();
            let output_paths = writer.write_all(&graph, &matrix, &gene_members)?;

            if self.config.enable_qc {
                let stats_path = self.config.output_dir.join("qc_stats.csv");
                if let Err(e) = crate::output::write_qc_stats(&qc_results, &stats_path) {
                    tracing::warn!("Failed to write QC stats: {}", e);
                }
            }

            tracing::info!("Prodigal pipeline complete. Output: {:?}", output_paths.output_dir);
            Ok(output_paths)
        }

        #[cfg(not(feature = "prodigal"))]
        {
            Err(Error::FeatureNotEnabled(
                "Prodigal gene calling requires the 'prodigal' feature flag. Rebuild with --features prodigal".into()
            ))
        }
    }

    /// Run pre-processing QC on all input files.
    fn run_qc(&self) -> Result<Vec<GenomeQC>> {
        let mut qc_results = Vec::new();

        // Try CheckM2
        let checkm_runner = self.config.checkm_database_path.clone()
            .map(|_| CheckmQcRunner::with_path("checkm2"))
            .or_else(|| {
                // Try to detect CheckM2
                CheckmQcRunner::detect()
            });

        // Run QC on each input file
        for input_path in &self.config.input_files {
            let genome_id = input_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();

            let mut qc = GenomeQC {
                genome_id,
                ..Default::default()
            };

            if let Some(runner) = &checkm_runner {
                let qc_runner: &dyn QcRunner = runner;
                if let Ok(qc_checkm) = qc_runner.run_qc(input_path) {
                    qc.completeness = qc_checkm.completeness;
                    qc.contamination = qc_checkm.contamination;
                    qc.genome_size = qc_checkm.genome_size;
                    qc.num_contigs = qc_checkm.num_contigs;
                    qc.n50 = qc_checkm.n50;
                    qc.mash_distance = qc_checkm.mash_distance;

                    // Update passed status based on CheckM metrics
                    let comp_threshold = self.config.qc_mode.min_completeness();
                    let cont_threshold = self.config.qc_mode.contamination_threshold();
                    qc.passed = qc.passed
                        && qc.completeness >= comp_threshold
                        && qc.contamination <= cont_threshold;
                }
            }

            qc_results.push(qc);
        }

        // Write QC output files
        let output_dir = &self.config.output_dir;
        fs::create_dir_all(output_dir)?;

        let stats_path = output_dir.join("qc_stats.csv");
        write_qc_stats(&qc_results, &stats_path)?;
        tracing::info!("Wrote QC statistics to {}", stats_path.display());

        let summary_path = output_dir.join("qc_summary.txt");
        write_qc_summary(&qc_results, &summary_path)?;
        tracing::info!("Wrote QC summary to {}", summary_path.display());

        Ok(qc_results)
    }

    /// Run GWAS analysis using pyseer.
    ///
    /// Generates input files from the pangenome graph and matrix,
    /// then invokes pyseer for gene-phenotype association testing.
    fn run_gwas(&self, graph: &PangenomeGraph, matrix: &BitPackedMatrix) -> Result<()> {
        let mut runner = PyseerRunner::new();

        // If a phenotype file is provided, use it
        if let Some(ref phenotype_path) = self.config.phenotype_file {
            if !phenotype_path.exists() {
                return Err(Error::InvalidInput(format!(
                    "Phenotype file not found: {}",
                    phenotype_path.display()
                )));
            }
            runner.with_phenotypes(phenotype_path.clone());
        }

        // If a distance matrix file exists from QC, use it
        let dist_path = self.config.output_dir.join("distance_matrix.csv");
        if dist_path.exists() {
            runner.with_distances(dist_path);
        }

        // Check if pyseer is available
        if !runner.is_available() {
            tracing::warn!("pyseer not installed. GWAS analysis skipped.");
            tracing::warn!("Install with: pip install pyseer");
            return Ok(());
        }

        // Run pyseer
        match runner.run(graph, matrix) {
            Ok(output) => {
                // Write GWAS results
                let gwas_path = self.config.output_dir.join("gwas_results.tsv");
                let mut file = std::fs::File::create(&gwas_path)?;
                use std::io::Write;
                writeln!(file, "cluster_id\teffect_size\tp_value\tfdr\tsignificant")?;
                for result in &output.results {
                    writeln!(
                        file,
                        "{}\t{:.6}\t{:.2e}\t{:.2e}\t{}",
                        result.snp_id,
                        result.effect_size,
                        result.p_value,
                        result.fdr,
                        if result.fdr < 0.05 { "yes" } else { "no" }
                    )?;
                }
                tracing::info!(
                    "GWAS: {} genes tested, {} significant (FDR < 0.05)",
                    output.snp_count,
                    output.significant_count
                );
                tracing::info!("Wrote GWAS results to {}", gwas_path.display());
            }
            Err(e) => {
                tracing::warn!("GWAS analysis failed: {}", e);
            }
        }

        Ok(())
    }

    /// Re-annotate input genomes with Bakta.
    ///
    /// Runs Bakta on FASTA/GenBank inputs and collects GFF3 output paths.
    /// GFF files are passed through unchanged.
    fn reannotate_inputs(&self, input_files: &[PathBuf]) -> Result<Vec<PathBuf>> {
        // Check if any files need re-annotation
        let needs_annotation = input_files.iter().any(|f: &PathBuf| {
            let ext = f.extension().and_then(|e: &std::ffi::OsStr| e.to_str()).unwrap_or("");
            !matches!(ext.to_lowercase().as_str(), "gff" | "gff3")
        });

        if !needs_annotation {
            tracing::info!("All input files are already in GFF3 format, skipping re-annotation");
            return Ok(input_files.to_vec());
        }

        // Try to detect Bakta
        let runner = match BaktaRunner::detect() {
            Some(runner) => runner,
            None => {
                tracing::warn!(
                    "Bakta not found. Install with: conda install -c conda-forge -c bioconda bakta"
                );
                tracing::warn!("Falling back to using input files directly");

                // Check if any files are GenBank format (can't be used without Bakta)
                for f in input_files {
                    if is_genbank_file(f) {
                        return Err(crate::Error::genbank_requires_bakta(f));
                    }
                }

                return Ok(input_files.to_vec());
            }
        };

        // Resolve database path
        let db_path = if let Some(ref path) = self.config.bakta_db_path {
            path.clone()
        } else {
            let resolved = BaktaRunner::resolve_db(None);
            if !resolved.exists() {
                if self.config.no_bakta_db_download {
                    return Err(crate::Error::bakta_db_not_found(&resolved));
                }
                // Auto-download the database
                let home = std::env::var("HOME")
                    .or_else(|_| std::env::var("USERPROFILE"))
                    .unwrap_or_else(|_| ".".to_string());
                let bakta_dir = PathBuf::from(home).join(".bakta");
                std::fs::create_dir_all(&bakta_dir)?;
                BaktaRunner::download_db(&bakta_dir, self.config.bakta_db_type)?
            } else {
                resolved
            }
        };

        let bakta_threads = if self.config.bakta_threads > 0 {
            self.config.bakta_threads
        } else {
            self.config.effective_threads()
        };

        // Use the detected bakta binary path and resolved db_path
        let annotator = BaktaRunner::new(runner.name_path(), db_path)
            .with_threads(bakta_threads)
            .with_output_dir(self.config.output_dir.clone())
            .with_keep_contig_headers(true);

        tracing::info!("Running Bakta re-annotation on {} files", input_files.len());

        let gff_paths = annotator.annotate_batch(input_files)?;

        tracing::info!("Re-annotation complete: {} GFF3 files", gff_paths.len());
        Ok(gff_paths)
    }

    /// Parse all input GFF3 files in parallel (using config input_files).
    ///
    /// Returns genes, genome IDs, and full contig DNA from the GFF FASTA section.
    #[allow(dead_code)]
    fn parse_inputs(&self) -> Result<(Vec<Gene>, Vec<GenomeId>, HashMap<(GenomeId, String), Vec<u8>>)> {
        self.parse_inputs_from(&self.config.input_files)
    }

    /// Parse input GFF3 files in parallel (from explicit file list).
    fn parse_inputs_from(&self, input_files: &[PathBuf]) -> Result<(Vec<Gene>, Vec<GenomeId>, HashMap<(GenomeId, String), Vec<u8>>)> {
        let results: Vec<(Vec<Gene>, GenomeId, HashMap<String, Vec<u8>>)> = input_files
            .par_iter()
            .map(|path: &PathBuf| {
                let genome_id = GenomeId::new(
                    path.file_stem()
                        .and_then(|s: &std::ffi::OsStr| s.to_str())
                        .unwrap_or("unknown")
                );

                match GffParser::open(path, genome_id.clone())
                    .and_then(|p: crate::io::GffParser| p.parse_genes_and_contigs())
                {
                    Ok((genes, contigs)) => {
                        if genes.is_empty() {
                            tracing::warn!("No genes found in {:?}", path);
                        }
                        (genes, genome_id, contigs)
                    }
                    Err(e) => {
                        tracing::error!("Failed to parse {:?}: {}. Genome will be excluded from analysis.", path, e);
                        (Vec::new(), genome_id, HashMap::new())
                    }
                }
            })
            .collect();

        let mut all_genes = Vec::new();
        let mut genome_ids = Vec::new();
        let mut all_contigs: HashMap<(GenomeId, String), Vec<u8>> = HashMap::new();

        for (genes, genome_id, contigs) in results {
            all_genes.extend(genes);
            genome_ids.push(genome_id.clone());
            for (contig_name, seq) in contigs {
                all_contigs.insert((genome_id.clone(), contig_name), seq);
            }
        }

        Ok((all_genes, genome_ids, all_contigs))
    }

    /// Cluster genes using MMseqs2-GPU or CPU fallback.
    fn cluster_genes(&self, genes: &[Gene]) -> Result<Vec<GeneCluster>> {
        // Try MMseqs2 first
        if self.config.enable_mmseqs && !self.config.force_cpu {
            if let Some(runner) = MMseqsRunner::detect() {
                tracing::info!("Using {} for clustering", runner.name());
                return runner.cluster(genes, self.config.cluster_identity, self.config.len_dif_percent);
            }
            tracing::info!("MMseqs2 not found, falling back to CPU clustering");
        }

        // CPU fallback
        let clusterer = CpuClusterer::new(self.config.effective_threads());
        tracing::info!("Using {} for clustering", clusterer.name());
        clusterer.cluster(genes, self.config.cluster_identity, self.config.len_dif_percent)
    }

    /// Build the pangenome graph from clusters, genes, and contig DNA.
    fn build_graph(&self, clusters: &[GeneCluster], genes: &[Gene], contig_dna: &HashMap<(GenomeId, String), Vec<u8>>) -> ConcurrentGraph {
        let builder = GraphBuilder::new()
            .with_min_support(self.config.min_support)
            .with_contig_dna(contig_dna);

        builder.build_concurrent(clusters, genes)
    }

    /// Run error correction on the graph.
    fn run_corrections(&self, graph: &ConcurrentGraph, num_genomes: usize) -> Result<()> {
        // Correction pipeline follows Panaroo's validated order:
        // 1. collapse_paralogs
        // 2. collapse_families (mistranslation, DNA, 0.98)
        // 3. collapse_families (families, protein, 0.7)
        // 4. trim_low_support_trailing_ends
        // 5. find_missing
        // 6. collapse_families (re-collapse)
        // 7. clean_misassembly_edges

        // Phase 4.0: Paralog resolution (must run before other corrections)
        tracing::info!("Phase 4.0: Running paralog resolution");
        let resolver = ParalogResolver::new();
        let paralog_stats = resolver.resolve(graph)?;
        tracing::info!(
            "Paralog resolution: {} paralogs detected, {} merges",
            paralog_stats.paralogs_detected,
            paralog_stats.nodes_merged
        );

        // Phase 4.1: Contamination removal (PanMiner extension — Panaroo lacks this step)
        let remover = ContaminationRemover::from_mode(&self.config.mode, num_genomes);
        remover.remove(graph)?;

        // Phase 4.2: Fragment merging with iterative multi-threshold collapsing
        // (matches Panaroo steps 2-3: collapse_families mistranslation + families)
        let merger = FragmentMerger::new()
            .with_collapse_thresholds(self.config.collapse_thresholds.clone());

        // Build sequences HashMap from graph nodes (all centroids)
        let sequences: std::collections::HashMap<String, Vec<u8>> = graph
            .nodes
            .iter()
            .flat_map(|entry| {
                let cluster_id = entry.key().clone();
                let centroids = entry.value().centroid_sequences.clone();
                centroids.into_iter().enumerate().filter_map(move |(i, seq)| {
                    if seq.is_empty() {
                        None
                    } else if i == 0 {
                        Some((cluster_id.to_string(), seq))
                    } else {
                        Some((format!("{}_centroid{}", cluster_id, i), seq))
                    }
                })
            })
            .collect();

        tracing::info!("Passing {} sequences to fragment merger", sequences.len());

        if sequences.is_empty() {
            tracing::warn!("No centroid sequences available for fragment merging");
        }

        // Create distance cache for reuse across correction passes (matches Panaroo Step 7->10)
        let mut distance_cache = DistanceCache::new();

        // Mistranslation correction at identity 0.99 (Panaroo step 2)
        merger.correct_mistranslations(graph, &sequences)?;

        // Rebuild sequences after mistranslation correction (graph was modified)
        let sequences: std::collections::HashMap<String, Vec<u8>> = graph
            .nodes
            .iter()
            .flat_map(|entry| {
                let cluster_id = entry.key().clone();
                let centroids = entry.value().centroid_sequences.clone();
                centroids.into_iter().enumerate().filter_map(move |(i, seq)| {
                    if seq.is_empty() {
                        None
                    } else if i == 0 {
                        Some((cluster_id.to_string(), seq))
                    } else {
                        Some((format!("{}_centroid{}", cluster_id, i), seq))
                    }
                })
            })
            .collect();

        // Iterative gene family collapsing from high to low threshold
        // (matches Panaroo step 3: collapse_families progressive collapsing)
        let mut total_collapsed = 0usize;
        for threshold in merger.collapse_thresholds() {
            let collapsed = merger.collapse_gene_families_with_threshold(
                graph, &sequences, *threshold, Some(&mut distance_cache)
            )?;
            total_collapsed += collapsed;
        }
        tracing::info!(
            "Iterative gene family collapsing: {} total merges across {} thresholds",
            total_collapsed,
            self.config.collapse_thresholds.len()
        );

        // Phase 4.3: Contig-end pruning (Panaroo step 4: trim_low_support_trailing_ends)
        // Must run AFTER collapsing so that merged nodes are properly evaluated
        let pruner = ContigEndPruner::from_mode(&self.config.mode, num_genomes);
        let pruning_stats = pruner.prune(graph)?;
        tracing::info!(
            "Contig-end pruning: removed {} nodes in {} iterations",
            pruning_stats.nodes_removed,
            pruning_stats.iterations
        );

        // Phase 4.4: Missing gene recovery (Panaroo step 5: find_missing)
        tracing::info!("Phase 4.4: Running missing gene recovery");
        self.run_missing_gene_recovery(graph)?;

        // Phase 4.5: Re-collapse families after missing gene recovery (Panaroo step 6)
        // Reuses distance cache to avoid recomputing known distances
        tracing::info!("Phase 4.5: Re-collapsing gene families after missing gene recovery");
        let sequences_after_recovery: std::collections::HashMap<String, Vec<u8>> = graph
            .nodes
            .iter()
            .flat_map(|entry| {
                let cluster_id = entry.key().clone();
                let centroids = entry.value().centroid_sequences.clone();
                centroids.into_iter().enumerate().filter_map(move |(i, seq)| {
                    if seq.is_empty() {
                        None
                    } else if i == 0 {
                        Some((cluster_id.to_string(), seq))
                    } else {
                        Some((format!("{}_centroid{}", cluster_id, i), seq))
                    }
                })
            })
            .collect();

        let merger2 = FragmentMerger::new()
            .with_collapse_thresholds(self.config.collapse_thresholds.clone());
        let mut total_recollapsed = 0usize;
        for threshold in merger2.collapse_thresholds() {
            let collapsed = merger2.collapse_gene_families_with_threshold(
                graph, &sequences_after_recovery, *threshold, Some(&mut distance_cache)
            )?;
            total_recollapsed += collapsed;
        }
        tracing::info!(
            "Post-recovery re-collapsing: {} total merges",
            total_recollapsed
        );

        // Phase 4.6: Misassembly edge cleaning (Panaroo step 7)
        let cleaner = MisassemblyEdgeCleaner::from_mode(&self.config.mode, num_genomes);
        let cleaning_stats = cleaner.clean(graph)?;
        tracing::info!("Misassembly edge cleaning: removed {} edges", cleaning_stats.edges_removed);

        Ok(())
    }

    /// Run missing gene recovery on the graph.
    ///
    /// This searches for genes that may have been missed during annotation
    /// by looking at flanking sequences around expected gene locations.
    ///
    /// For each edge in the graph, checks if any genome is missing
    /// one of the connected genes. If so, searches the contig sequences
    /// for a match using semi-global alignment.
    fn run_missing_gene_recovery(&self, graph: &ConcurrentGraph) -> Result<()> {
        // Extract cluster sequences from graph nodes
        let cluster_sequences: HashMap<String, Vec<u8>> = graph
            .nodes
            .iter()
            .filter_map(|entry| {
                let cluster_id = entry.key().as_str().to_string();
                let node = entry.value();
                node.centroid_sequences.first().map(|seq| (cluster_id.clone(), seq.clone()))
            })
            .collect();

        if cluster_sequences.is_empty() {
            tracing::info!("Missing gene recovery: no cluster sequences available");
            return Ok(());
        }

        // Extract contig sequences organized by genome.
        // Each node stores contig_sequences as HashMap<contig_name, Sequence>.
        // We need to map these to the genomes that have this node so that
        // when searching for a missing gene in a specific genome, we search
        // the correct contigs.
        let mut genome_contig_sequences: HashMap<GenomeId, Vec<Vec<u8>>> = HashMap::new();
        for entry in graph.nodes.iter() {
            let node = entry.value();
            for genome_id in &node.genomes {
                for (_contig_name, seq) in &node.contig_sequences {
                    genome_contig_sequences
                        .entry(genome_id.clone())
                        .or_default()
                        .push(seq.clone());
                }
            }
        }

        // Also build a flat map for the recoverer (it just needs contig sequences to search)
        // Key by genome so we only search relevant contigs
        let contig_sequences: HashMap<String, Vec<u8>> = genome_contig_sequences
            .iter()
            .flat_map(|(genome_id, seqs)| {
                seqs.iter().enumerate().map(move |(i, seq)| {
                    (format!("{}_contig{}", genome_id.as_str(), i), seq.clone())
                })
            })
            .collect();

        if contig_sequences.is_empty() {
            tracing::info!("Missing gene recovery: no contig sequences available");
            return Ok(());
        }

        tracing::info!(
            "Missing gene recovery: {} cluster sequences, {} contig sequences from {} genomes",
            cluster_sequences.len(),
            contig_sequences.len(),
            genome_contig_sequences.len()
        );

        // Perform missing gene recovery
        let remove_by_consensus = matches!(self.config.mode, CorrectionMode::Strict);
        let recoverer = MissingGeneRecoverer::new()
            .with_min_identity(0.70)
            .with_search_window(5000)
            .with_remove_by_consensus(remove_by_consensus);

        let stats = recoverer.recover(graph, &contig_sequences, &cluster_sequences)?;
        tracing::info!(
            "Missing gene recovery: {} genes recovered, {} nodes removed by consensus",
            stats.genes_recovered,
            stats.nodes_removed_by_consensus
        );

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
            .with_chunk_size(2) // 3 files, chunk size 2 -> 2 chunks
            .with_enable_qc(false); // Disable QC for test with mock files

        let pipeline = PanminerPipeline::new(config);
        let paths = pipeline.run().expect("Pipeline should run successfully with chunks");

        // The matrix should have 3 genomes
        let matrix_content = std::fs::read_to_string(paths.matrix_csv.as_ref().unwrap()).unwrap();
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
