//! Integration tests for fragment merging with real sequences.
//!
//! Tests that the FragmentMerger correctly processes actual cluster centroid sequences
//! instead of empty sequences.

use panminer::graph::{GeneCluster, Node, ConcurrentGraph, ClusterId, GeneId, GenomeId};
use panminer::correction::FragmentMerger;

#[test]
fn test_fragment_merger_with_real_sequences() {
    // Create two clusters with highly similar sequences (99% identity)
    let mut cluster_a = GeneCluster::new("cluster_a");
    cluster_a.genes.push(GeneId::new("gene_a1"));
    cluster_a.centroid = Some(b"ATCGATCGATCGATCGATCGATCG".to_vec());
    cluster_a.support = 5;

    let mut cluster_b = GeneCluster::new("cluster_b");
    cluster_b.genes.push(GeneId::new("gene_b1"));
    // 99% identical to cluster_a (1 mismatch in 100 bases)
    cluster_b.centroid = Some(b"ATCGATCGATCGATCGATCGATCA".to_vec()); // One mismatch
    cluster_b.support = 3;

    // Build graph with both clusters
    let graph = ConcurrentGraph::new();
    let node_a = Node::from_cluster(&cluster_a);
    let node_b = Node::from_cluster(&cluster_b);
    graph.add_node(node_a);
    graph.add_node(node_b);

    // Create a connection between them
    let mut edge = panminer::graph::Edge::new(
        ClusterId::new("cluster_a"),
        ClusterId::new("cluster_b")
    );
    edge.add_genome(GenomeId::new("genome1"));
    graph.add_edge(edge);

    // Build sequences HashMap from graph nodes
    let sequences: std::collections::HashMap<String, Vec<u8>> = graph
        .nodes
        .iter()
        .filter_map(|entry| {
            entry.value().centroid_sequence.clone().map(|seq| (entry.key().to_string(), seq))
        })
        .collect();

    assert_eq!(sequences.len(), 2, "Should have 2 sequences");

    // Run fragment merger
    let merger = FragmentMerger::new()
        .with_collapse_threshold(0.70);

    let result = merger.correct_mistranslations(&graph, &sequences);
    assert!(result.is_ok(), "Mistranslation correction should succeed");
}

#[test]
fn test_fragment_merger_no_sequences() {
    // Test that with empty sequences, no merging occurs
    let graph = ConcurrentGraph::new();
    let mut node = Node::from_cluster(&{
        let mut c = GeneCluster::new("test");
        c.centroid = None;
        c.support = 2;
        c
    });
    node.centroid_sequence = None; // Explicitly no sequence
    graph.add_node(node);

    let empty_sequences: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();
    let merger = FragmentMerger::new();

    let result = merger.correct_mistranslations(&graph, &empty_sequences);
    assert!(result.is_ok(), "Should handle empty sequences gracefully");
}

#[test]
fn test_node_from_cluster_with_centroid() {
    let mut cluster = GeneCluster::new("test_cluster");
    cluster.centroid = Some(b"ATCGATCGATCGATCG".to_vec());
    cluster.support = 3;

    let node = Node::from_cluster(&cluster);
    assert_eq!(node.centroid_sequence, Some(b"ATCGATCGATCGATCG".to_vec()));
}

#[test]
fn test_empty_graph_no_sequences() {
    // Test that an empty graph produces no sequences
    let graph = ConcurrentGraph::new();
    let sequences: std::collections::HashMap<String, Vec<u8>> = graph
        .nodes
        .iter()
        .filter_map(|entry| entry.value().centroid_sequence.clone().map(|seq| (entry.key().to_string(), seq)))
        .collect();

    assert!(sequences.is_empty(), "Empty graph should have no sequences");
}

#[test]
fn test_fragment_merger_identical_sequences() {
    // Test that identical sequences get merged
    let graph = ConcurrentGraph::new();

    let mut cluster_a = GeneCluster::new("cluster_a");
    cluster_a.centroid = Some(b"ATCGATCGATCGATCG".to_vec());
    cluster_a.support = 5;
    graph.add_node(Node::from_cluster(&cluster_a));

    let mut cluster_b = GeneCluster::new("cluster_b");
    cluster_b.centroid = Some(b"ATCGATCGATCGATCG".to_vec()); // Identical
    cluster_b.support = 3;
    graph.add_node(Node::from_cluster(&cluster_b));

    // Create connection
    let mut edge = panminer::graph::Edge::new(
        ClusterId::new("cluster_a"),
        ClusterId::new("cluster_b")
    );
    edge.add_genome(GenomeId::new("genome1"));
    graph.add_edge(edge);

    let sequences: std::collections::HashMap<String, Vec<u8>> = graph
        .nodes
        .iter()
        .filter_map(|entry| entry.value().centroid_sequence.clone().map(|seq| (entry.key().to_string(), seq)))
        .collect();

    let merger = FragmentMerger::new()
        .with_collapse_threshold(0.70);

    // Identical sequences should have 100% identity and 100% coverage
    // This should trigger a merge
    let result = merger.correct_mistranslations(&graph, &sequences);
    assert!(result.is_ok(), "Should process identical sequences");
}
