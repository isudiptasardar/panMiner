---
name: current_state
description: PanMiner current state and features
type: project
---

**Current Status (2026-04-10):**
- Rust v0.1.0 with cargo check passing (minimal warnings for optional feature code)
- 100 tests passing: 81 unit + 5 fragment merge + 9 integration + 4 P1 integration + 5 doc

**Implemented P1 Features (Complete):**
- Contig-end pruning (`src/correction/contig_end.rs`) - Iterative removal of degree-1 nodes at contig ends
- Structural variant matrix (`src/output/sv_matrix.rs`) - Gene triplet TSV output
- Expanded integration tests (`tests/p1_integration_tests.rs`) - 4 new integration tests

**Implemented P2 Features (Complete):**
- GWAS Integration (`src/gwas/pyseer.rs`) - Pyseer wrapper with distance/phenotype file generation
- Real MSA Output (`src/output/alignment.rs`) - MAFFT/Clustal/PRANK subprocess invocation
- Parquet Output (`src/output/parquet.rs`) - Columnar format with optional feature flag
- HTML Visualization (`src/output/html_viz.rs`) - Interactive d3.js force-directed graph

**Previously Implemented Features (Complete):**
- GFF3/FASTA parsing with memory-mapped I/O
- MMseqs2 clustering (GPU detection, CPU fallback)
- Graph construction with DashMap (concurrent, lock-free)
- Contamination removal (iterative degree-1 node removal)
- Fragment merging with real centroid sequences
  - Mistranslation correction (95% coverage, 99% identity)
  - Gene family collapsing (70% identity, shared neighbor detection)
- Missing gene recovery (k-mer search, 5000bp flanking, 70% threshold)
- MSA integration (MAFFT, Clustal Omega, PRANK)
- Output formats: CSV, TSV, GML, JSON, JSONL, QC stats
- Streaming pipeline with chunked processing and zstd compression
- CheckM2 pre-processing QC

**Configuration:**
- Correction modes: strict (≥5%), default (≥1%), sensitive (no deletion)
- Identity threshold: configurable (default 0.98)
- Clustering: MMseqs2-GPU (auto-detected) or CPU fallback
- Output formats: matrix, alignment, graph, json, parquet, html

**Known Warnings (Expected):**
- `HtmlVizWriter` shows dead_code warning when `viz` feature disabled (expected behavior)

**Data Structures:**
- `Node.contig_sequences`: HashMap<String, Sequence> - stores contig sequences for each node
- `Node.is_contig_end`: bool - marks if node is at contig end
- `GeneCluster.centroid`: Option<Sequence> - stores the representative sequence
- `PangenomeGraph.nodes`: HashMap<ClusterId, Node> - standard graph representation
- `ConcurrentGraph.nodes`: DashMap<ClusterId, Node> - concurrent graph representation

**Error Correction Pipeline:**
1. Contamination removal - Removes low-support degree-1 nodes
2. Contig-end pruning - Removes degree-1 nodes at contig ends (P1)
3. Fragment merging - Merges genes with ≥95% coverage and ≥99% identity
4. Gene family collapsing - Merges gene families at ≥70% identity
5. Missing gene recovery - K-mer search (11-mer, 70% threshold) in 5000bp flanking regions
