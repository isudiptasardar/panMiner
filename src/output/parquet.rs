//! Parquet output format support for PanMiner.
//!
//! This module provides Parquet file output for the pangenome graph and matrix,
//! enabling efficient columnar storage and integration with data analytics tools.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{BooleanArray, Int32Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;

use crate::error::{Error, Result};
use crate::graph::{PangenomeGraph, BitPackedMatrix};

/// Writer for Parquet output format.
pub struct ParquetWriter;

impl ParquetWriter {
    /// Create a new Parquet writer.
    pub fn new() -> Self {
        ParquetWriter
    }

    /// Write the presence/absence matrix to a Parquet file.
    ///
    /// The matrix is stored as a table where:
    /// - Each row represents a cluster
    /// - Each column represents a genome (binary presence/absence)
    ///
    /// # Arguments
    ///
    /// * `matrix` - The bit-packed presence/absence matrix
    /// * `path` - Output file path
    pub fn write_matrix(&self, matrix: &BitPackedMatrix, path: &Path) -> Result<()> {
        let num_genomes = matrix.num_genomes();
        let num_clusters = matrix.num_clusters();

        // Build schema for the matrix
        let mut fields: Vec<Field> = vec![Field::new("cluster_id", DataType::Utf8, false)];

        for genome_name in &matrix.genome_names {
            fields.push(Field::new(genome_name.as_str(), DataType::Boolean, false));
        }

        let schema = Arc::new(Schema::new(fields));

        // Build record batch
        let mut column_data: Vec<Arc<dyn arrow::array::Array>> = Vec::new();

        // Cluster IDs column
        let cluster_ids: Vec<&str> = matrix.cluster_ids.iter().map(|s| s.as_str()).collect();
        column_data.push(Arc::new(StringArray::from(cluster_ids)));

        // For each genome, create a boolean array
        for genome_idx in 0..num_genomes {
            let values: Vec<bool> = (0..num_clusters)
                .map(|cluster_idx| matrix.get(genome_idx, cluster_idx))
                .collect();
            column_data.push(Arc::new(BooleanArray::from(values)));
        }

        let batch = RecordBatch::try_new(schema.clone(), column_data)
            .map_err(|e| Error::Arrow(e.to_string()))?;

        // Write to file
        let file = std::fs::File::create(path)?;
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))
            .map_err(|e| Error::Parquet(e.to_string()))?;
        writer.write(&batch).map_err(|e| Error::Parquet(e.to_string()))?;
        writer.close().map_err(|e| Error::Parquet(e.to_string()))?;

        Ok(())
    }

    /// Write the pangenome graph to Parquet format.
    ///
    /// This creates two Parquet files:
    /// - `nodes.parquet`: Node data (cluster_id, support, genomes, etc.)
    /// - `edges.parquet`: Edge data (from, to, genomes, support)
    ///
    /// # Arguments
    ///
    /// * `graph` - The pangenome graph
    /// * `path` - Output file path
    pub fn write_graph(&self, graph: &PangenomeGraph, path: &Path) -> Result<()> {
        // Write nodes
        self.write_nodes(graph, path)?;

        // Write edges
        self.write_edges(graph, path)?;

        Ok(())
    }

    /// Write nodes to Parquet format.
    fn write_nodes(&self, graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let schema = Arc::new(Schema::new(vec![
            Field::new("cluster_id", DataType::Utf8, false),
            Field::new("support", DataType::Int32, false),
            Field::new("is_paralog", DataType::Boolean, false),
            Field::new("genome_count", DataType::Int32, false),
        ]));

        let mut cluster_ids: Vec<&str> = Vec::new();
        let mut supports: Vec<i32> = Vec::new();
        let mut is_paralogs: Vec<bool> = Vec::new();
        let mut genome_counts: Vec<i32> = Vec::new();

        for (cluster_id, node) in &graph.nodes {
            cluster_ids.push(cluster_id.as_str());
            supports.push(node.support as i32);
            is_paralogs.push(node.is_paralog);
            genome_counts.push(node.genomes.len() as i32);
        }

        let column_data: Vec<Arc<dyn arrow::array::Array>> = vec![
            Arc::new(StringArray::from(cluster_ids)),
            Arc::new(Int32Array::from(supports)),
            Arc::new(BooleanArray::from(is_paralogs)),
            Arc::new(Int32Array::from(genome_counts)),
        ];

        let batch = RecordBatch::try_new(schema, column_data)
            .map_err(|e| Error::Arrow(e.to_string()))?;

        // Write to nodes.parquet
        let node_path = path.parent()
            .map(|p| p.join("nodes.parquet"))
            .unwrap_or_else(|| PathBuf::from("nodes.parquet"));

        let file = std::fs::File::create(&node_path)?;
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))
            .map_err(|e| Error::Parquet(e.to_string()))?;
        writer.write(&batch).map_err(|e| Error::Parquet(e.to_string()))?;
        writer.close().map_err(|e| Error::Parquet(e.to_string()))?;

        Ok(())
    }

    /// Write edges to Parquet format.
    fn write_edges(&self, graph: &PangenomeGraph, path: &Path) -> Result<()> {
        let schema = Arc::new(Schema::new(vec![
            Field::new("from", DataType::Utf8, false),
            Field::new("to", DataType::Utf8, false),
            Field::new("support", DataType::Int32, false),
            Field::new("genome_count", DataType::Int32, false),
        ]));

        let mut from_ids: Vec<&str> = Vec::new();
        let mut to_ids: Vec<&str> = Vec::new();
        let mut supports: Vec<i32> = Vec::new();
        let mut genome_counts: Vec<i32> = Vec::new();

        for ((from, to), edge) in &graph.edges {
            from_ids.push(from.as_str());
            to_ids.push(to.as_str());
            supports.push(edge.support as i32);
            genome_counts.push(edge.genomes.len() as i32);
        }

        let column_data: Vec<Arc<dyn arrow::array::Array>> = vec![
            Arc::new(StringArray::from(from_ids)),
            Arc::new(StringArray::from(to_ids)),
            Arc::new(Int32Array::from(supports)),
            Arc::new(Int32Array::from(genome_counts)),
        ];

        let batch = RecordBatch::try_new(schema, column_data)
            .map_err(|e| Error::Arrow(e.to_string()))?;

        // Write to edges.parquet
        let edge_path = path.parent()
            .map(|p| p.join("edges.parquet"))
            .unwrap_or_else(|| PathBuf::from("edges.parquet"));

        let file = std::fs::File::create(&edge_path)?;
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))
            .map_err(|e| Error::Parquet(e.to_string()))?;
        writer.write(&batch).map_err(|e| Error::Parquet(e.to_string()))?;
        writer.close().map_err(|e| Error::Parquet(e.to_string()))?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Node, GeneCluster, GenomeId, ClusterId};
    use std::fs::File;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parquet_matrix_output() {
        let mut matrix = BitPackedMatrix::new(3, 2);
        matrix.set_genome_names(vec!["g1".to_string(), "g2".to_string(), "g3".to_string()]);
        matrix.set_cluster_ids(vec!["c1".to_string(), "c2".to_string()]);

        matrix.set(0, 0, true);
        matrix.set(1, 0, true);
        matrix.set(2, 1, true);

        let mut temp = NamedTempFile::new().unwrap();
        let writer = ParquetWriter::new();

        let result = writer.write_matrix(&matrix, temp.path());
        if let Err(err) = &result {
            eprintln!("Matrix output failed: {}", err);
        }
        result.unwrap();
    }

    #[test]
    fn test_parquet_graph_output() {
        let mut graph = PangenomeGraph::new();

        let node1 = Node::from_cluster(&{
            let mut c = GeneCluster::new("c1");
            c.support = 5;
            c
        });
        graph.add_node(node1);

        let node2 = Node::from_cluster(&{
            let mut c = GeneCluster::new("c2");
            c.support = 3;
            c
        });
        graph.add_node(node2);

        let mut edge = crate::graph::Edge::new(ClusterId::new("c1"), ClusterId::new("c2"));
        edge.add_genome(GenomeId::new("genome1"));
        graph.add_edge(edge);

        let temp = NamedTempFile::new().unwrap();
        let writer = ParquetWriter::new();

        let result = writer.write_graph(&graph, temp.path());
        if let Err(err) = &result {
            eprintln!("Graph output failed: {}", err);
        }
        result.unwrap();
    }
}
