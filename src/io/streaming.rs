//! Streaming pipeline for bounded memory usage with large datasets.

use std::path::{Path, PathBuf};
use std::collections::HashMap;
use rayon::prelude::*;

use crate::config::PanminerConfig;
use crate::error::{Error, Result};
use crate::graph::{Gene, GeneCluster, GenomeId};
use super::gff::GffParser;
use super::compress::{write_compressed, read_compressed};

/// Partial graph built from a chunk of genomes.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PartialGraph {
    /// Chunk ID
    pub chunk_id: usize,
    /// Gene adjacencies (edges) as ClusterIds
    pub adjacencies: Vec<(String, String, String)>, // (contig, from_cluster, to_cluster)
    /// Genome IDs in this chunk
    pub genome_ids: Vec<String>,
}

/// Streaming pipeline for processing large datasets.
pub struct StreamingPipeline {
    chunk_size: usize,
    intermediate_dir: PathBuf,
}

impl StreamingPipeline {
    pub fn new(config: PanminerConfig) -> Self {
        let chunk_size = config.chunk_size;
        let intermediate_dir = config.temp_dir.join("panminer_intermediate");
        Self { chunk_size, intermediate_dir }
    }

    pub fn process_chunks_with_clusters(
        &self,
        genome_files: &[PathBuf],
        clusters: &[GeneCluster],
    ) -> Result<Vec<PathBuf>> {
        std::fs::create_dir_all(&self.intermediate_dir)?;

        let gene_to_cluster: HashMap<String, String> = clusters
            .iter()
            .flat_map(|c| c.genes.iter().map(move |g| (g.to_string(), c.id.to_string())))
            .collect();

        genome_files
            .chunks(self.chunk_size)
            .enumerate()
            .map(|(chunk_id, chunk)| {
                self.process_chunk(chunk_id, chunk, &gene_to_cluster)
            })
            .collect()
    }

    fn process_chunk(
        &self,
        chunk_id: usize,
        genome_files: &[PathBuf],
        gene_to_cluster: &HashMap<String, String>,
    ) -> Result<PathBuf> {
        let all_genes: Vec<Gene> = genome_files
            .par_iter()
            .flat_map(|path| {
                let genome_id = GenomeId::new(path.file_stem().and_then(|s| s.to_str()).unwrap_or("unknown"));
                GffParser::open(path, genome_id).and_then(|p| p.parse_genes()).unwrap_or_default()
            })
            .collect();

        let mut contig_genes: HashMap<(String, String), Vec<&Gene>> = HashMap::new();
        for gene in &all_genes {
            let key = (gene.genome_id.to_string(), gene.contig.clone());
            contig_genes.entry(key).or_default().push(gene);
        }

        let mut adjacencies = Vec::new();
        for ((_genome, contig), mut genes) in contig_genes {
            genes.sort_by_key(|g| g.start);
            for window in genes.windows(2) {
                if let (Some(c1), Some(c2)) = (
                    gene_to_cluster.get(&window[0].id.to_string()),
                    gene_to_cluster.get(&window[1].id.to_string()),
                ) {
                    adjacencies.push((contig.clone(), c1.clone(), c2.clone()));
                }
            }
        }

        let partial = PartialGraph {
            chunk_id,
            adjacencies,
            genome_ids: genome_files.iter().map(|p| p.file_stem().unwrap().to_str().unwrap().to_string()).collect(),
        };

        let path = self.intermediate_dir.join(format!("chunk_{}.zst", chunk_id));
        let serialized = bincode::serialize(&partial).map_err(|e| Error::Output(e.to_string()))?;
        write_compressed(&path, &serialized)?;
        Ok(path)
    }

    pub fn load_intermediate(&self, path: &Path) -> Result<PartialGraph> {
        let compressed = read_compressed(path)?;
        bincode::deserialize(&compressed).map_err(|e| Error::Output(e.to_string()))
    }
}
