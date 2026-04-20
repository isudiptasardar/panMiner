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
    /// (genome_id, from_cluster, to_cluster)
    pub adjacencies: Vec<(String, String, String)>,
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
                match GffParser::open(path, genome_id) {
                    Ok(parser) => match parser.parse_genes() {
                        Ok(genes) => genes,
                        Err(e) => {
                            tracing::warn!("Failed to parse genes from {:?}: {}", path, e);
                            Vec::new()
                        }
                    },
                    Err(e) => {
                        tracing::warn!("Failed to open GFF file {:?}: {}", path, e);
                        Vec::new()
                    }
                }
            })
            .collect();

        let mut contig_genes: HashMap<(String, String), Vec<&Gene>> = HashMap::new();
        for gene in &all_genes {
            let key = (gene.genome_id.to_string(), gene.contig.clone());
            contig_genes.entry(key).or_default().push(gene);
        }

        let mut adjacencies = Vec::new();
        for ((genome, _contig), mut genes) in contig_genes {
            genes.sort_by_key(|g| g.start);
            for window in genes.windows(2) {
                if let (Some(c1), Some(c2)) = (
                    gene_to_cluster.get(&window[0].id.to_string()),
                    gene_to_cluster.get(&window[1].id.to_string()),
                ) {
                    adjacencies.push((genome.clone(), c1.clone(), c2.clone()));
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_streaming_pipeline_creation() {
        let config = PanminerConfig::new()
            .with_input_files(vec![PathBuf::from("test.gff")])
            .with_output_dir(PathBuf::from("test_output"));
        let pipeline = StreamingPipeline::new(config);
        assert_eq!(pipeline.chunk_size, 100); // default chunk_size
    }

    #[test]
    fn test_streaming_pipeline_custom_chunk_size() {
        let config = PanminerConfig::new()
            .with_input_files(vec![PathBuf::from("test.gff")])
            .with_output_dir(PathBuf::from("test_output"))
            .with_chunk_size(50);
        let pipeline = StreamingPipeline::new(config);
        assert_eq!(pipeline.chunk_size, 50);
    }

    #[test]
    fn test_partial_graph_serialization() {
        let partial = PartialGraph {
            chunk_id: 0,
            adjacencies: vec![
                ("contig1".to_string(), "cluster_A".to_string(), "cluster_B".to_string()),
                ("contig1".to_string(), "cluster_B".to_string(), "cluster_C".to_string()),
            ],
            genome_ids: vec!["genome1".to_string(), "genome2".to_string()],
        };

        let serialized = bincode::serialize(&partial).unwrap();
        let deserialized: PartialGraph = bincode::deserialize(&serialized).unwrap();
        assert_eq!(deserialized.chunk_id, 0);
        assert_eq!(deserialized.adjacencies.len(), 2);
        assert_eq!(deserialized.genome_ids.len(), 2);
    }

    #[test]
    fn test_partial_graph_empty() {
        let partial = PartialGraph {
            chunk_id: 1,
            adjacencies: vec![],
            genome_ids: vec![],
        };

        let serialized = bincode::serialize(&partial).unwrap();
        let deserialized: PartialGraph = bincode::deserialize(&serialized).unwrap();
        assert_eq!(deserialized.chunk_id, 1);
        assert!(deserialized.adjacencies.is_empty());
        assert!(deserialized.genome_ids.is_empty());
    }

    #[test]
    fn test_partial_graph_round_trip_compressed() {
        let dir = tempfile::tempdir().unwrap();
        let partial = PartialGraph {
            chunk_id: 42,
            adjacencies: vec![
                ("contig1".to_string(), "A".to_string(), "B".to_string()),
            ],
            genome_ids: vec!["g1".to_string()],
        };

        let path = dir.path().join("chunk_42.zst");
        let serialized = bincode::serialize(&partial).unwrap();
        write_compressed(&path, &serialized).unwrap();

        let compressed = read_compressed(&path).unwrap();
        let deserialized: PartialGraph = bincode::deserialize(&compressed).unwrap();
        assert_eq!(deserialized.chunk_id, 42);
        assert_eq!(deserialized.adjacencies.len(), 1);
        assert_eq!(deserialized.genome_ids[0], "g1");
    }
}
