# PanMiner P1 Features Implementation Design

> **Date**: 2026-04-10  
> **Author**: Brainstorming session with user input  
> **Status**: Approved for implementation

---

## Overview

This document describes the implementation plan for PanMiner P1 (High Priority) features to achieve Panaroo feature parity:

1. **Contig-end pruning** - Recursive degree-1 node removal at contig ends
2. **Structural variant matrix** - Gene triplet presence/absence output
3. **Expanded integration tests** - End-to-end pipeline testing with real data scenarios

---

## 1. Contig-End Pruning

### 1.1 Feature Description

Panaroo has a specific "contig-end pruning" step that removes fragmented genes at contig boundaries. This is distinct from general contamination removal (which removes low-support degree-1 nodes).

**Current PanMiner behavior:**
- Contamination remover removes low-support degree-1 nodes
- No specific handling for contig-end fragmentation

**Required behavior:**
- Detect nodes that are at contig ends (single gene on contig)
- Iteratively remove these degree-1 nodes
- Configurable threshold for minimum support

### 1.2 Algorithm

```
function contig_end_pruning(graph, min_support=1):
    changed = true
    while changed:
        changed = false
        for node in graph.nodes:
            if node.is_at_contig_end() and node.support < min_support:
                graph.remove_node(node)
                changed = true
```

### 1.3 Data Structure Changes

Add field to `Node`:
```rust
pub struct Node {
    // ... existing fields ...
    pub is_contig_end: bool,  // Marked during graph building
}
```

### 1.4 New Module: `correction/contig_end.rs`

```rust
pub struct ContigEndPruner {
    min_support: usize,
}

impl ContigEndPruner {
    pub fn new(min_support: usize) -> Self;
    pub fn prune(&self, graph: &ConcurrentGraph) -> Result<PruningStats>;
}
```

---

## 2. Structural Variant Matrix

### 2.1 Feature Description

Panaroo outputs a "structural variant matrix" - gene triplet presence/absence across genomes. This captures co-occurrence patterns of nearby genes.

**Current:** None (only presence/absence matrix)

**Required:** Gene triplet matrix in TSV format

### 2.2 Algorithm

For each node in the graph:
1. Find neighboring nodes (via edges)
2. For each pair of neighbors, check which genomes contain both
3. Output as gene triplets: `clusterA_clusterB_genome1, clusterA_clusterB_genome2, ...`

### 2.3 Output Format

```tsv
## Structural Variant Matrix (gene triplets)
## Format: clusterA_clusterB = presence/absence across genomes
clusterA_clusterB	genome1	genome2	genome3
gene1_gene2	1	1	0
gene1_gene3	1	0	1
```

### 2.4 New Module: `output/structural_variants.rs`

```rust
pub struct StructuralVariantMatrix {
    triplets: Vec<(ClusterId, ClusterId, Vec<bool>)>,
}

impl StructuralVariantMatrix {
    pub fn from_graph(graph: &PangenomeGraph) -> Self;
    pub fn write_tsv(&self, path: &Path) -> Result<()>;
}
```

---

## 3. Expanded Integration Tests

### 3.1 Test Coverage

| Test Name | Description |
|-----------|-------------|
| `test_pipeline_contig_end_pruning` | Verifies degree-1 node removal at contig ends |
| `test_pipeline_sv_matrix` | Verifies structural variant matrix output |
| `test_pipeline_with_paralogs` | Tests paralog handling |
| `test_pipeline_large_dataset` | Tests with 10+ genomes |

### 3.2 Test Data

Use real GFF3 files or generate synthetic data with:
- Multiple genes per contig
- Known adjacency patterns
- Control paralog scenarios

---

## Implementation Plan

### Phase 1: Contig-End Pruning
1. Add `is_contig_end` field to `Node` struct
2. Mark contig ends during `GraphBuilder::build_concurrent`
3. Create `ContigEndPruner` with iterative removal
4. Wire into pipeline after contamination removal

### Phase 2: Structural Variant Matrix
1. Implement `StructuralVariantMatrix` struct
2. Algorithm to extract gene triplets from graph
3. TSV writer with proper headers
4. Integrate into `OutputWriter`

### Phase 3: Expanded Tests
1. Create test data generator for contig-end scenarios
2. Test contig-end pruning with known graph topologies
3. Test structural variant matrix output
4. Large-scale integration test

---

## Error Handling

All new modules should use `thiserror` for typed errors:

```rust
#[derive(thiserror::Error, Debug)]
pub enum ContigEndError {
    #[error("Graph is empty")]
    EmptyGraph,
    #[error("Node not found: {0}")]
    NodeNotFound(ClusterId),
}
```

---

## Testing Strategy

- Unit tests for each new module
- Integration tests for end-to-end scenarios
- Compare output with Panaroo on same test data (when possible)

---

## Success Criteria

- All P1 features implemented and passing tests
- Performance within 2x of Panaroo on comparable datasets
- Output format compatible with Panaroo where specified
