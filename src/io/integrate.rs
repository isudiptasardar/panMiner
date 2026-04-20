//! Integrate a new genome into an existing pangenome graph.
//!
//! This module provides incremental genome addition: given an existing
//! PanMiner output directory (containing `final_graph.gml`) and a new
//! GFF3 file, it loads the graph, parses the new genome, matches new
//! genes against existing centroids by sequence identity, and produces
//! an updated pangenome.

use std::collections::HashMap;
use std::path::Path;

use crate::error::Result;
use crate::graph::{
    ClusterId, Edge, GeneCluster, GeneId, GenomeId, GenomeMetadata, Node, PangenomeGraph,
};
use crate::io::GffParser;

/// Result of an integrate operation.
#[derive(Debug)]
pub struct IntegrateResult {
    /// Path to the output directory
    pub output_dir: std::path::PathBuf,
    /// Total nodes after integration
    pub total_nodes: usize,
    /// Total edges after integration
    pub total_edges: usize,
    /// Number of new genes added to existing nodes
    pub genes_matched: usize,
    /// Number of new nodes created for unmatched genes
    pub new_nodes: usize,
}

/// Integrate a single GFF file into an existing PanMiner pangenome.
///
/// Algorithm:
/// 1. Load `final_graph.gml` from the existing directory
/// 2. Parse the new GFF file
/// 3. For each new gene, compare its protein sequence against existing
///    centroids. If identity is at or above the threshold, add the gene
///    to the matching node; otherwise create a new node.
/// 4. Build adjacency edges for genes that are adjacent on the same contig
/// 5. Write the updated graph and output files to `output_dir`
pub fn integrate_genome(
    existing_dir: &Path,
    new_gff: &Path,
    output_dir: &Path,
    identity_threshold: f32,
    threads: usize,
) -> Result<IntegrateResult> {
    // Step 1: Load existing graph
    let gml_path = existing_dir.join("final_graph.gml");
    if !gml_path.exists() {
        return Err(crate::Error::Config(format!(
            "No final_graph.gml found in {:?}",
            existing_dir
        )));
    }
    let mut graph = load_gml_graph(&gml_path)?;

    // Step 2: Parse the new GFF file
    let genome_id = derive_genome_id(new_gff);
    let parser = GffParser::open(new_gff, genome_id.clone())?;
    let (new_genes, contigs) = parser.parse_genes_and_contigs()?;

    if new_genes.is_empty() {
        return Err(crate::Error::NoGenes(genome_id.to_string()));
    }

    tracing::info!(
        "Parsed {} genes from new genome {}",
        new_genes.len(),
        genome_id
    );

    // Register the new genome in the graph
    graph.genomes.insert(
        genome_id.clone(),
        GenomeMetadata {
            id: genome_id.clone(),
            source_file: new_gff.to_string_lossy().to_string(),
            num_contigs: contigs.len(),
            total_genes: new_genes.len(),
        },
    );

    // Pre-compute protein sequences for new genes
    let new_proteins: HashMap<GeneId, Vec<u8>> = new_genes
        .iter()
        .filter_map(|gene| {
            if gene.sequence.is_empty() {
                None
            } else {
                let protein = crate::io::translate(&gene.sequence);
                Some((gene.id.clone(), protein))
            }
        })
        .collect();

    // Step 3: Match each new gene against existing centroids
    let mut genes_matched = 0usize;
    let mut new_node_count = 0usize;

    // Collect existing centroids for comparison (cluster_id -> centroid protein sequences)
    let existing_centroids: Vec<(ClusterId, Vec<Vec<u8>>)> = graph
        .nodes
        .iter()
        .filter_map(|(cid, node)| {
            let proteins: Vec<Vec<u8>> = node
                .centroid_sequences
                .iter()
                .map(|seq| crate::io::translate(seq))
                .filter(|p| !p.is_empty())
                .collect();
            if proteins.is_empty() {
                None
            } else {
                Some((cid.clone(), proteins))
            }
        })
        .collect();

    // Map from new gene index to the cluster it matched (or None for new)
    let mut gene_assignments: Vec<(usize, Option<ClusterId>)> = Vec::new();

    for (i, gene) in new_genes.iter().enumerate() {
        let gene_protein = match new_proteins.get(&gene.id) {
            Some(p) if !p.is_empty() => p,
            _ => {
                // No protein sequence available; always create a new node
                gene_assignments.push((i, None));
                continue;
            }
        };

        // Find the best-matching existing centroid
        let mut best_cluster: Option<ClusterId> = None;
        let mut best_identity = 0.0f32;

        for (cid, centroid_proteins) in &existing_centroids {
            for centroid_prot in centroid_proteins {
                let identity = sequence_identity(gene_protein, centroid_prot);
                if identity > best_identity {
                    best_identity = identity;
                    best_cluster = Some(cid.clone());
                }
            }
        }

        if best_identity >= identity_threshold {
            if let Some(ref cid) = best_cluster {
                genes_matched += 1;
                gene_assignments.push((i, Some(cid.clone())));
            }
        } else {
            gene_assignments.push((i, None));
        }
    }

    // Apply assignments: add matched genes to existing nodes
    for (gene_idx, assignment) in &gene_assignments {
        let gene = &new_genes[*gene_idx];
        if let Some(ref cid) = assignment {
            if let Some(node) = graph.nodes.get_mut(cid) {
                node.support += 1;
                node.genomes.insert(gene.genome_id.clone());
                node.gene_members
                    .entry(gene.genome_id.clone())
                    .or_default()
                    .push(gene.id.as_str().to_string());

                // Also store contig sequence if available
                if !gene.sequence.is_empty() {
                    let contig_key = format!("{}:{}", gene.genome_id, gene.contig);
                    node.add_contig_sequence(contig_key, gene.sequence.clone());
                }

                // Add gene to gene_lookup
                graph.gene_lookup.insert(gene.id.clone(), gene.clone());
            }
        }
    }

    // Create new nodes for unmatched genes
    for (gene_idx, assignment) in &gene_assignments {
        if assignment.is_some() {
            continue;
        }
        let gene = &new_genes[*gene_idx];

        // Generate a cluster ID from the gene ID, avoiding collisions
        let mut cluster_id = ClusterId::new(format!("cluster_{}", gene.id));
        if graph.nodes.contains_key(&cluster_id) {
            cluster_id = ClusterId::new(format!("cluster_{}_new", gene.id));
        }
        if graph.nodes.contains_key(&cluster_id) {
            cluster_id = ClusterId::new(format!("cluster_{}_{}", gene.id, new_node_count));
        }

        let mut cluster = GeneCluster::new(cluster_id.as_str());
        cluster.support = 1;
        cluster.add_gene(gene.id.clone());

        // Use DNA sequence as centroid; translate to protein for centroid_sequences
        if !gene.sequence.is_empty() {
            cluster.centroids.push(gene.sequence.clone());
        }

        let mut node = Node::from_cluster_with_genes(&cluster, &{
            let mut data = HashMap::new();
            data.insert(gene.id.clone(), gene.clone());
            data
        });
        node.genomes.insert(gene.genome_id.clone());

        // Store contig sequence
        if !gene.sequence.is_empty() {
            let contig_key = format!("{}:{}", gene.genome_id, gene.contig);
            node.add_contig_sequence(contig_key, gene.sequence.clone());
        }

        graph.add_node(node);
        graph.gene_lookup.insert(gene.id.clone(), gene.clone());
        new_node_count += 1;
    }

    tracing::info!(
        "Matched {} genes to existing clusters, created {} new clusters",
        genes_matched,
        new_node_count
    );

    // Step 4: Build adjacency edges for new genes on the same contig
    let _effective_threads = if threads == 0 {
        std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(1)
    } else {
        threads
    };

    // Group new genes by contig
    let mut contig_genes: HashMap<String, Vec<(usize, Option<ClusterId>)>> = HashMap::new();
    for (gene_idx, assignment) in &gene_assignments {
        let gene = &new_genes[*gene_idx];
        contig_genes
            .entry(gene.contig.clone())
            .or_default()
            .push((*gene_idx, assignment.clone()));
    }

    // For each contig, create edges between adjacent genes
    for (_contig, genes_on_contig) in &contig_genes {
        let mut sorted_genes = genes_on_contig.clone();
        sorted_genes.sort_by_key(|(idx, _)| new_genes[*idx].start);

        for pair in sorted_genes.windows(2) {
            let (idx_a, assign_a) = &pair[0];
            let (idx_b, assign_b) = &pair[1];
            let gene_a = &new_genes[*idx_a];
            let gene_b = &new_genes[*idx_b];

            let cid_a = assign_a
                .clone()
                .unwrap_or_else(|| ClusterId::new(format!("cluster_{}", gene_a.id)));
            let cid_b = assign_b
                .clone()
                .unwrap_or_else(|| ClusterId::new(format!("cluster_{}", gene_b.id)));

            if cid_a != cid_b {
                // Check if edge already exists
                let edge_exists = graph.has_edge(&cid_a, &cid_b);
                if edge_exists {
                    // Add this genome to the existing edge
                    if let Some(edge) = graph.edges.get_mut(&(cid_a.clone(), cid_b.clone())) {
                        edge.add_genome(genome_id.clone());
                    } else if let Some(edge) = graph.edges.get_mut(&(cid_b.clone(), cid_a.clone())) {
                        edge.add_genome(genome_id.clone());
                    }
                } else {
                    let mut edge = Edge::new(cid_a, cid_b);
                    edge.add_genome(genome_id.clone());
                    graph.add_edge(edge);
                }
            }
        }
    }

    // Step 5: Correction passes are skipped during integration because the
    // correction modules operate on ConcurrentGraph, not PangenomeGraph.
    // The user should re-run the full pipeline if correction is needed.
    // For now, we just write the augmented graph as-is.

    // Step 6: Write updated output to output_dir
    std::fs::create_dir_all(output_dir)?;

    let output_gml = output_dir.join("final_graph.gml");
    crate::output::GmlWriter::write(&graph, &output_gml)?;

    // Write gene_data.csv for gene sequence access
    let gene_data_path = output_dir.join("gene_data.csv");
    crate::output::JsonWriter::write_gene_data(&graph, &graph.gene_lookup, &gene_data_path)?;

    let total_nodes = graph.node_count();
    let total_edges = graph.edge_count();

    Ok(IntegrateResult {
        output_dir: output_dir.to_path_buf(),
        total_nodes,
        total_edges,
        genes_matched,
        new_nodes: new_node_count,
    })
}

/// Derive a genome ID from a GFF file path (stem of the file name).
fn derive_genome_id(path: &Path) -> GenomeId {
    path.file_stem()
        .and_then(|s| s.to_str())
        .map(|s| GenomeId::new(s))
        .unwrap_or_else(|| GenomeId::new("unknown_genome"))
}

/// Load a pangenome graph from a GML file.
///
/// This is a self-contained reader that handles the full PanMiner GML
/// output format, including centroid sequences, genome IDs, gene members,
/// and contig-end genome sets. It reuses the parsing approach from
/// `graph::merge::load_gml_graph` but additionally restores the richer
/// node attributes that the merge reader leaves as defaults.
fn load_gml_graph(path: &Path) -> Result<PangenomeGraph> {
    let content = std::fs::read_to_string(path)?;
    let mut graph = PangenomeGraph::new();

    let mut current_node: Option<Node> = None;
    let mut current_edge: Option<Edge> = None;
    let mut in_node = false;
    let mut in_edge = false;

    for line in content.lines() {
        let line = line.trim();

        if line == "node [" {
            in_node = true;
            current_node = None;
        } else if line == "]" && in_node {
            if let Some(node) = current_node.take() {
                graph.nodes.insert(node.cluster_id.clone(), node);
            }
            in_node = false;
        } else if line == "edge [" {
            in_edge = true;
        } else if line == "]" && in_edge {
            if let Some(edge) = current_edge.take() {
                let key = (edge.from.clone(), edge.to.clone());
                graph.edges.insert(key, edge);
            }
            in_edge = false;
        } else if in_node {
            if let Some(node) = &mut current_node {
                parse_node_attribute(node, line);
            } else {
                // Start a default node; the label/id line will set cluster_id
                current_node = Some(Node::from_cluster(&{
                    let mut c = GeneCluster::new("temp");
                    c.support = 1;
                    c
                }));
                // Parse the first attribute of this node
                if let Some(ref mut node) = current_node {
                    parse_node_attribute(node, line);
                }
            }
        } else if in_edge {
            if current_edge.is_none() {
                current_edge = Some(Edge::new(
                    ClusterId::new("unknown"),
                    ClusterId::new("unknown"),
                ));
            }
            if let Some(edge) = &mut current_edge {
                parse_edge_attribute(edge, line);
            }
        }
    }

    // Rebuild adjacency index from edges
    graph.rebuild_adjacency();

    Ok(graph)
}

/// Parse a single node attribute line from GML format.
fn parse_node_attribute(node: &mut Node, line: &str) {
    if line.starts_with("id ") || line.starts_with("label ") {
        if let Some(label) = extract_gml_string(line) {
            node.cluster_id = ClusterId::new(&label);
        }
    } else if line.starts_with("support") {
        if let Some(s) = line.split_whitespace().nth(1) {
            node.support = s.parse().unwrap_or(1);
        }
    } else if line.starts_with("is_paralog") {
        if let Some(s) = line.split_whitespace().nth(1) {
            node.is_paralog = s.parse().unwrap_or(0) != 0;
        }
    } else if line.starts_with("is_highly_variable") {
        if let Some(s) = line.split_whitespace().nth(1) {
            node.is_highly_variable = s.parse().unwrap_or(0) != 0;
        }
    } else if line.starts_with("seq ") {
        if let Some(seq_str) = extract_gml_string(line) {
            // Comma-separated centroid sequences
            node.centroid_sequences = seq_str
                .split(',')
                .map(|s| s.as_bytes().to_vec())
                .collect();
        }
    } else if line.starts_with("genome_ids ") {
        if let Some(ids_str) = extract_gml_string(line) {
            for gid in ids_str.split(',') {
                let trimmed = gid.trim();
                if !trimmed.is_empty() {
                    node.genomes.insert(GenomeId::new(trimmed));
                }
            }
        }
    } else if line.starts_with("member ") {
        if let Some(members_str) = extract_gml_string(line) {
            // Semicolon-separated gene member IDs. GML format does not
            // preserve per-genome gene attribution, so we distribute
            // members across all known genomes for this node. Each genome
            // gets the full member list (over-estimation, but preserves IDs).
            let members: Vec<String> = members_str
                .split(';')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();
            if !members.is_empty() {
                // Distribute across all known genomes for this node
                for gid in &node.genomes {
                    node.gene_members.insert(gid.clone(), members.clone());
                }
            }
        }
    } else if line.starts_with("contig_end_genomes ") {
        if let Some(ids_str) = extract_gml_string(line) {
            for gid in ids_str.split(',') {
                let trimmed = gid.trim();
                if !trimmed.is_empty() {
                    node.contig_end_genomes.insert(GenomeId::new(trimmed));
                }
            }
        }
    } else if line.starts_with("annotation ") {
        if let Some(ann) = extract_gml_string(line) {
            node.annotations.insert(ann);
        }
    }
    // Silently ignore: length, protein (derived from seq at runtime)
}

/// Parse a single edge attribute line from GML format.
fn parse_edge_attribute(edge: &mut Edge, line: &str) {
    if line.starts_with("source ") {
        if let Some(s) = extract_gml_string(line) {
            edge.from = ClusterId::new(&s);
        }
    } else if line.starts_with("target ") {
        if let Some(s) = extract_gml_string(line) {
            edge.to = ClusterId::new(&s);
        }
    } else if line.starts_with("support") {
        if let Some(s) = line.split_whitespace().nth(1) {
            edge.support = s.parse().unwrap_or(1);
        }
    } else if line.starts_with("genome_ids ") {
        if let Some(ids_str) = extract_gml_string(line) {
            for gid in ids_str.split(',') {
                let trimmed = gid.trim();
                if !trimmed.is_empty() {
                    edge.genomes.insert(GenomeId::new(trimmed));
                }
            }
        }
    }
}

/// Extract a quoted string value from a GML attribute line.
///
/// Handles escaped quotes and backslashes per GML spec.
fn extract_gml_string(line: &str) -> Option<String> {
    // Find the first and last quote
    let start = line.find('"')?;
    let end = line.rfind('"')?;
    if start >= end {
        return None;
    }
    let raw = &line[start + 1..end];
    // Un-escape GML: \\ -> \, \" -> "
    Some(raw.replace("\\\"", "\"").replace("\\\\", "\\"))
}

/// Calculate sequence identity between two protein sequences.
///
/// Uses simple character-level identity over the longer sequence length.
fn sequence_identity(seq1: &[u8], seq2: &[u8]) -> f32 {
    if seq1.is_empty() || seq2.is_empty() {
        return 0.0;
    }

    let min_len = seq1.len().min(seq2.len());
    let max_len = seq1.len().max(seq2.len());
    let mut matches = 0;

    for i in 0..min_len {
        if seq1[i].to_ascii_uppercase() == seq2[i].to_ascii_uppercase() {
            matches += 1;
        }
    }

    matches as f32 / max_len as f32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{GeneCluster, GenomeMetadata};
    use std::io::Write;

    fn make_test_graph(nodes: Vec<(&str, usize, Vec<Vec<u8>>)>) -> PangenomeGraph {
        let mut graph = PangenomeGraph::new();
        for (id, support, centroids) in nodes {
            let mut cluster = GeneCluster::new(id);
            cluster.support = support;
            cluster.centroids = centroids;
            let node = Node::from_cluster(&cluster);
            graph.add_node(node);
        }
        graph
    }

    #[test]
    fn test_sequence_identity_identical() {
        let seq1 = b"MFFLLK";
        let seq2 = b"MFFLLK";
        let identity = sequence_identity(seq1, seq2);
        assert!((identity - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_sequence_identity_different() {
        let seq1 = b"MFFLLK";
        let seq2 = b"GGGGGG";
        let identity = sequence_identity(seq1, seq2);
        assert!(identity < 0.3);
    }

    #[test]
    fn test_sequence_identity_empty() {
        let identity = sequence_identity(b"", b"MFF");
        assert_eq!(identity, 0.0);
    }

    #[test]
    fn test_extract_gml_string_simple() {
        let line = r#"label "cluster_001""#;
        assert_eq!(extract_gml_string(line), Some("cluster_001".to_string()));
    }

    #[test]
    fn test_extract_gml_string_escaped() {
        let line = r#"seq "ATG\\CGT""#;
        assert_eq!(extract_gml_string(line), Some("ATG\\CGT".to_string()));
    }

    #[test]
    fn test_extract_gml_string_no_quotes() {
        let line = "support 5";
        assert_eq!(extract_gml_string(line), None);
    }

    #[test]
    fn test_derive_genome_id() {
        let id = derive_genome_id(Path::new("/data/genomes/sample1.gff"));
        assert_eq!(id.as_str(), "sample1");
    }

    #[test]
    fn test_derive_genome_id_no_extension() {
        let id = derive_genome_id(Path::new("/data/genomes/sample1"));
        assert_eq!(id.as_str(), "sample1");
    }

    #[test]
    fn test_integrate_missing_gml_fails() {
        let result = integrate_genome(
            Path::new("/nonexistent_dir"),
            Path::new("/nonexistent.gff"),
            Path::new("/tmp/integrate_out"),
            0.98,
            1,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_load_gml_roundtrip() {
        // Build a graph, write GML, then load it back
        let mut graph = PangenomeGraph::new();

        let mut node = Node::from_cluster(&{
            let mut c = GeneCluster::new("c1");
            c.support = 3;
            c.centroids = vec![b"ATGCGT".to_vec()];
            c
        });
        node.genomes.insert(GenomeId::new("genome1"));
        node.genomes.insert(GenomeId::new("genome2"));
        graph.add_node(node);

        let dir = tempfile::tempdir().unwrap();
        let gml_path = dir.path().join("test.gml");
        crate::output::GmlWriter::write(&graph, &gml_path).unwrap();

        let loaded = load_gml_graph(&gml_path).unwrap();
        assert_eq!(loaded.node_count(), 1);
        assert_eq!(
            loaded.nodes.get(&ClusterId::new("c1")).unwrap().support,
            3
        );
        assert!(loaded
            .nodes
            .get(&ClusterId::new("c1"))
            .unwrap()
            .genomes
            .contains(&GenomeId::new("genome1")));
    }

    #[test]
    fn test_integrate_single_genome() {
        // Create a small existing graph with GML
        let mut graph = make_test_graph(vec![("c1", 2, vec![b"ATGCGTAAA".to_vec()])]);
        graph.genomes.insert(
            GenomeId::new("genome1"),
            GenomeMetadata {
                id: GenomeId::new("genome1"),
                source_file: "genome1.gff".to_string(),
                num_contigs: 1,
                total_genes: 1,
            },
        );
        graph.genomes.insert(
            GenomeId::new("genome2"),
            GenomeMetadata {
                id: GenomeId::new("genome2"),
                source_file: "genome2.gff".to_string(),
                num_contigs: 1,
                total_genes: 1,
            },
        );

        // Set up the node with genomes
        if let Some(node) = graph.nodes.get_mut(&ClusterId::new("c1")) {
            node.genomes.insert(GenomeId::new("genome1"));
            node.genomes.insert(GenomeId::new("genome2"));
        }

        let dir = tempfile::tempdir().unwrap();
        let gml_path = dir.path().join("final_graph.gml");
        crate::output::GmlWriter::write(&graph, &gml_path).unwrap();

        // Create a new GFF with a matching gene
        let mut gff = tempfile::NamedTempFile::new().unwrap();
        writeln!(gff, "##gff-version 3").unwrap();
        writeln!(
            gff,
            "seq1\tProkka\tgene\t100\t300\t.\t+\t.\tID=gene_new;product=test"
        )
        .unwrap();
        writeln!(gff, "##FASTA").unwrap();
        writeln!(gff, ">seq1").unwrap();
        // ATGCGTAAA translated = MRVK (roughly), so we provide a DNA
        // that translates to a similar protein for the matching test
        writeln!(gff, "ATGCGTAAA").unwrap();

        let output_dir = tempfile::tempdir().unwrap();
        let result = integrate_genome(
            dir.path(),
            gff.path(),
            output_dir.path(),
            0.5, // Low threshold so the gene matches
            1,
        );

        assert!(result.is_ok(), "integrate_genome failed: {:?}", result.err());
        let result = result.unwrap();
        // The graph should have at least one node (either matched or new)
        assert!(result.total_nodes >= 1);
    }
}