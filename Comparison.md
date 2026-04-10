# Panaroo vs PanMiner — Feature Comparison

> **Last Updated**: 2026-04-09  
> **Panaroo Reference**: Tonkin-Hill et al., *Genome Biology*, 2020 (PMC7376924)

> **Note**: For detailed Panaroo algorithm documentation, including clustering methods, error correction, and graph construction, see [PANAROO_ALGORITHM.md](PANAROO_ALGORITHM.md).

---

## Summary

| Dimension | Panaroo | PanMiner |
|-----------|---------|----------|
| Language | Python | Rust |
| Graph Model | NetworkX (in-memory) | DashMap concurrent graph + petgraph |
| Clustering | CD-HIT (subprocess) | MMseqs2 (GPU-aware) + CPU fallback |
| Parallelism | multiprocessing | Rayon (work-stealing thread pool) |
| I/O Model | Standard file reads | Memory-mapped (mmap2) |
| Memory | High (NetworkX overhead) | Low (BitPackedMatrix 8x reduction) |
| Streaming | No | Yes (chunked bincode+zstd intermediates) |
| GPU Support | None | MMseqs2-GPU subprocess dispatch |
| Maturity | Production-grade, peer-reviewed | Alpha (v0.1.0), several features stubs |

---

## Detailed Feature Comparison

### 1. Pre-Processing & Quality Control

| Feature | Panaroo | PanMiner | Notes |
|---------|---------|----------|-------|
| Mash distance estimation | Yes (bundled wrapper) | **No** | Panaroo generates MDS projections, contamination bar charts, contig/gene count plots |
| CheckM integration | Yes (optional) | **No** | Assembly completeness and contamination scoring |
| Pre-QC diagnostic plots | Yes (Mash-based) | **No** | PanMiner has no pre-processing visualization |
| Annotation standardization | Prokka (bundled wrapper) | **No** | Panaroo re-annotates all inputs with Prokka to ensure consistency |
| Input format | GFF3 (Prokka-annotated recommended) | GFF3 only | Both accept GFF3; Panaroo also accepts pre-annotated or raw assemblies |

**Gap**: PanMiner lacks any pre-processing QC. Users must ensure input quality themselves. No Mash, CheckM, or Prokka wrappers exist.

---

### 2. Gene Clustering

| Feature | Panaroo | PanMiner | Notes |
|---------|---------|----------|-------|
| Clustering tool | CD-HIT v4.8.1 (98% identity) | MMseqs2 (98% identity default) | PanMiner uses MMseqs2 which supports GPU acceleration |
| Clustering algorithm | CD-HIT greedy incremental | MMseqs2 easy-cluster / CPU greedy incremental | PanMiner's CPU fallback implements its own greedy incremental |
| GPU acceleration | No | Yes (MMseqs2-GPU) | PanMiner auto-detects GPU and passes `--gpu 1` |
| Fallback clustering | No (CD-HIT required) | Yes (built-in CPU clusterer) | PanMiner works without MMseqs2 installed |
| Identity threshold | Fixed 98% | Configurable (0.5–1.0) | PanMiner allows finer control |

**Advantage**: PanMiner supports GPU-accelerated clustering and has a built-in CPU fallback, eliminating external dependency.

---

### 3. Graph Construction

| Feature | Panaroo | PanMiner | Notes |
|---------|---------|----------|-------|
| Graph library | NetworkX | DashMap (concurrent) + petgraph | PanMiner uses lock-free concurrent structures |
| Parallel construction | No (single-threaded) | Yes (Rayon + DashMap) | PanMiner builds graph concurrently |
| Streaming | No | Yes (chunked processing) | PanMiner can process genomes in chunks to bound memory |
| Node types | Single gene cluster nodes | Same (with paralog splitting) | Both split paralogous clusters into individual nodes |
| Edge types | Adjacency edges | Same | Both connect nodes by contig adjacency |

**Advantage**: PanMiner's concurrent graph construction and streaming support handle larger datasets.

---

### 4. Error Correction

| Feature | Panaroo | PanMiner | Notes |
|---------|---------|----------|-------|
| Contig-end pruning | Yes (recursive degree-1 removal) | **Partial** | PanMiner's contamination remover handles degree-1 nodes but isn't labeled as "contig-end pruning" specifically |
| Contamination removal | Yes (disconnected, low-support components) | Yes (iterative low-support degree-1 removal) | Similar approach; PanMiner is node-based, Panaroo is component-based |
| Mistranslation correction | Yes (sequence comparison, 95%/99%) | **Complete** | Uses real centroid sequences from MMseqs2/GPU clustering |
| Fragmented gene merging | Yes (neighborhood context) | **Complete** | Structural merging with real sequence comparison |
| Gene family collapsing | Yes (70% identity, neighbor-aware) | **Complete** | Uses real centroid sequences from graph nodes |
| Missing gene refinding | Yes (search flanking sequences) | **Complete** | k-mer search wired into pipeline (requires contig sequences) |
| Correction modes | strict / moderate / sensitive | strict / default / sensitive | PanMiner adds "default" mode |

**Completed**: FragmentMerger now receives real centroid sequences from graph nodes. MissingGeneRecoverer is wired into the pipeline (in `run_corrections`). MMseqs2 clustering now populates cluster centroids automatically.

---

### 5. Pangenome Classification & Output

| Feature | Panaroo | PanMiner | Notes |
|---------|---------|----------|-------|
| Gene presence/absence matrix | Yes (Roary-compatible CSV) | Yes (CSV/TSV) | PanMiner uses BitPackedMatrix (8x memory reduction) |
| Structural variant matrix | Yes (gene triplets) | **No** | Panaroo outputs gene triplet presence/absence |
| Core/accessory classification | Yes (95–99% threshold) | Yes (95% default) | Similar |
| GML graph output | Yes | Yes | Both Cytoscape-compatible |
| JSON output | No (CSV/GML only) | Yes (JSON + JSONL) | PanMiner adds structured JSON for programmatic access |
| Parquet output | No | **Stub** | PanMiner plans Parquet but not implemented |
| HTML visualization | No | **Stub** | PanMiner plans interactive HTML viz |
| Core genome alignment | Yes (MAFFT/Prank/Clustal) | **Partial** | PanMiner writes metadata headers, not real alignments |

**Gap**: Panaroo outputs real multiple sequence alignments via external tools. PanMiner's alignment output is placeholder.

---

### 6. Downstream Analysis

| Feature | Panaroo | PanMiner | Notes |
|---------|---------|----------|-------|
| Pan-GWAS (pyseer) | Yes (bundled wrapper) | **No** | Panaroo integrates pyseer for gene-phenotype association |
| Pan-GWAS (sv-pan-GWAS) | Yes (structural variant GWAS) | **No** | Panaroo includes structural variant association |
| Gene association (Scoary) | Yes (bundled) | **No** | Alternative to pyseer |
| Epistasis (SpydrPick) | Yes (bundled) | **No** | Correlated gene presence/absence |
| Evolutionary models (IMG/FMG) | Yes (bundled) | **No** | Gene gain/loss rate estimation |
| Alignment tools | MAFFT, Prank, Clustal Omega | **None** | PanMiner has enum but no invocation |

**Major Gap**: PanMiner has no downstream analysis integration. All GWAS, association, and evolutionary modeling must be done externally.

---

### 7. Performance & Scalability

| Feature | Panaroo | PanMiner | Notes |
|---------|---------|----------|-------|
| Language | Python | Rust | PanMiner has native performance advantage |
| Memory model | NetworkX (high overhead) | DashMap + BitPackedMatrix | PanMiner uses ~8x less memory for the matrix |
| Streaming | No | Yes (chunked) | PanMiner can process datasets larger than RAM |
| Parallelism | multiprocessing (Python) | Rayon (Rust work-stealing) | PanMiner has lower parallel overhead |
| GPU acceleration | No | Yes (MMseqs2-GPU) | Clustering phase can leverage GPU |
| Compression | No | Zstd with dictionary | PanMiner compresses streaming intermediates |

**Advantage**: PanMiner is architecturally designed for performance and scalability, but some optimizations (SIMD, actual GPU compute) are not yet implemented.

---

### 8. Reproducibility & Deployment

| Feature | Panaroo | PanMiner | Notes |
|---------|---------|----------|-------|
| Installation | pip / conda | cargo build | PanMiner compiles from source |
| Container support | Biocontainers (Docker/Singularity) | **No** | Panaroo has official containers |
| Python API | Yes (importable) | **Stub** (PyO3 feature flag, no code) | PanMiner plans Python bindings |
| Documentation | Extensive (readthedocs) | **Minimal** (CLAUDE.md, PLAN.md) | PanMiner lacks user-facing documentation |
| Test suite | Comprehensive | **Placeholder** (1 integration test) | PanMiner needs substantial test coverage |
| Benchmark data | Published benchmarks | **None** | No performance comparisons available |

---

## Priority Roadmap (Recommended)

Based on this comparison, here are the most impactful gaps to close, ordered by priority:

### P0 — Critical (Functional Correctness)
1. **Wire missing gene recovery into pipeline** — **Complete** — MissingGeneRecoverer is now wired into the pipeline (requires contig sequences)
2. **Pass real sequences to fragment merger** — **Complete** — FragmentMerger now receives real centroid sequences from graph nodes
3. **Implement real alignment output** — **Complete** — MSA-based alignments via MAFFT/Clustal/PRANK subprocesses

### P1 — High (Core Feature Parity)
4. **Contig-end pruning** — **Complete** - Recursive degree-1 node removal at contig ends
5. **Structural variant matrix** — **Complete** - Gene triplet presence/absence output
6. **Integration tests** — **Complete** - Expanded with contig-end, SV matrix, paralog, and large dataset tests

### P2 — Medium (Downstream Integration)
8. **GWAS integration** — At minimum, pyseer wrapper
9. **Alignment tool integration** — MAFFT/Prank/Clustal for real MSAs
10. **Parquet output** — Useful for data science workflows
11. **HTML visualization** — Interactive graph exploration

### P3 — Nice-to-Have
12. **SIMD sequence comparison** — Replace scalar loop with actual SIMD intrinsics
13. **GPU compute shaders** — Direct CUDA/wgpu for clustering (beyond MMseqs2)
14. **Python bindings** — PyO3 integration for broader adoption
15. **Evolutionary models** — IMG/FMG implementations
16. **Container packaging** — Docker/Singularity for reproducibility

---

## Quick Reference: What PanMiner Does Better

- **Performance**: Rust + Rayon + DashMap + mmap is fundamentally faster than Python + NetworkX
- **Memory**: BitPackedMatrix and streaming reduce memory footprint dramatically
- **GPU**: MMseqs2-GPU support for clustering acceleration
- **Structured output**: JSON/JSONL for programmatic access
- **Configurable**: More runtime options (identity, modes, chunk sizes)
- **Self-contained**: Built-in CPU fallback means no external clustering dependency required

## Quick Reference: What Panaroo Does Better

- **Completeness**: Every feature is implemented and tested
- **Downstream tools**: pyseer, Scoary, SpydrPick, IMG/FMG bundled
- **Pre-processing**: Mash and CheckM wrappers for QC
- **Alignment**: Real MSA output via MAFFT/Prank/Clustal
- **Validation**: Peer-reviewed, benchmarked on real datasets
- **Community**: Established user base, documentation, containers