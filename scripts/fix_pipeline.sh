#!/bin/bash
cat << 'INNER_EOF' > pipeline_patch.rs
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
INNER_EOF

# Extract lines 1-54 of pipeline.rs
head -n 54 src/pipeline.rs > pipeline_new.rs
# Append patch
cat pipeline_patch.rs >> pipeline_new.rs
# Extract from "tracing::info!("Phase 3: Building pangenome graph (in-memory)");" downwards
sed -n '/tracing::info!("Phase 3: Building pangenome graph (in-memory)");/,$p' src/pipeline.rs >> pipeline_new.rs
mv pipeline_new.rs src/pipeline.rs
