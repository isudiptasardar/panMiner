# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**IMPORTANT: Before developing, editing, or writing any code, ALWAYS read these files first:**
- **[README.md](README.md)** — Installation and usage instructions
- **[Specs.md](Specs.md)** — Complete feature specifications, configuration reference, and feature status matrix
- **[Comparison.md](Comparison.md)** — Panaroo vs PanMiner feature comparison with priority roadmap
- **[Architecture.md](Architecture.md)** — Folder structure, module logic, data flow, and design decisions
- **[PANAROO_ALGORITHM.md](PANAROO_ALGORITHM.md)** — Panaroo algorithm reference for PanMiner development

These documents define the current state of the project, known gaps, and planned direction. All development decisions must account for the information in these files.

## Project Overview

PanMiner (panminer) is a modern pangenome analysis tool written in Rust. It processes GFF3 files to build pangenome graphs, with support for GPU-accelerated clustering via MMseqs2, pre-processing QC via CheckM2, and CPU fallback. The project follows a layered architecture inspired by Panaroo but with modern optimizations (memory-mapped I/O, DashMap concurrent graphs, Rayon parallelism).

## Quick Start

### Installation

```bash
# From crates.io
cargo install panminer

# From source
cargo install --path .
```

See [README.md](README.md) for detailed installation options (conda, manual, etc.).

**GPU Detection:** The installation script (`bash install.sh`) automatically detects NVIDIA GPUs and offers to install MMseqs2 with GPU support. CheckM2 detection is automatic for pre-processing QC.

**QC Options:** Use `--no-qc` to disable pre-processing QC, `--qc-mode strict|default|sensitive` to set stringency, or `--checkm-database` to specify CheckM2 database path.

### Build and Test Commands

```bash
# Build in dev mode
cargo build

# Build optimized release
cargo build --release

# Run all tests
cargo test

# Run specific test
cargo test test_name

# Run tests with all features
cargo test --features full

# Check without building
cargo check

# Display CLI help
cargo run -- --help

# Run with inputs
cargo run -- input1.gff input2.gff -o output_dir
```

## Feature Flags

- `cpu` (default): CPU processing support
- `mmseqs`: MMseqs2 integration for clustering
- `parquet`: Parquet output format support
- `python`: PyO3 Python bindings
- `viz`: HTML visualization output
- `full`: All features enabled

## Current Status

| Category | Status |
|----------|--------|
| Compilation | ✅ Compiles with no errors |
| Tests | ✅ 134 tests passing (126 unit + 8 Bakta integration) |
| Documentation | ✅ README, Specs, Comparison, Architecture, CLAUDE.md complete |
| GPU Detection | ✅ Automatic detection during installation |
| Critical Gaps | See [Comparison.md](Comparison.md#priority-roadmap) |

### Known Critical Gaps

These are the most important functional gaps that must be addressed (see [Comparison.md](Comparison.md) for full details):

1. **Mash distance estimation** — Panaroo generates MDS projections using Mash, not yet implemented
2. **No downstream analysis** — No GWAS, alignment tools, or evolutionary models

## Architecture

```
src/
├── lib.rs              # Library entry point, re-exports
├── main.rs             # CLI with clap
├── config.rs           # Configuration structs (CorrectionMode, OutputFormat)
├── error.rs            # Error types (thiserror-based)
├── pipeline.rs         # Main pipeline orchestration (PanminerPipeline)
│
├── io/                 # Layer 1: I/O & Memory
│   ├── mmap.rs         # Memory-mapped file wrapper
│   ├── gff.rs          # GFF3 parser (mmap-based)
│   ├── fasta.rs        # FASTA parser
│   ├── compress.rs     # Zstd compression utilities
│   ├── streaming.rs    # Chunked processing for large datasets
│   ├── bakta.rs        # Bakta annotation runner (subprocess)
│   └── qc_traits.rs    # QC runner traits (CheckM2)
│
├── clustering/         # Layer 2: Compute
│   ├── traits.rs       # Clusterer trait
│   ├── mmseqs.rs       # MMseqs2-GPU integration
│   └── cpu.rs          # CPU fallback with SIMD
│
├── graph/              # Layer 3: Data Structures
│   ├── types.rs        # Gene, GeneCluster, Node, Edge, PangenomeGraph
│   ├── concurrent.rs   # DashMap-based graph (ConcurrentGraph)
│   ├── matrix.rs       # BitPackedMatrix (8x memory reduction)
│   └── builder.rs      # Graph construction
│
├── correction/         # Layer 2: Error Correction
│   ├── contamination.rs
│   ├── fragment.rs
│   ├── missing.rs
│   ├── paralog.rs
│   └── misassembly.rs
│
└── output/             # Layer 4: Output Generation
    ├── matrix.rs       # CSV/TSV presence-absence
    ├── alignment.rs    # FASTA alignments
    ├── graph.rs        # GML output
    └── json.rs         # JSON/JSONL output
```

See [Architecture.md](Architecture.md) for full module-by-module logic, data flow diagrams, and design decisions.

## Pipeline Flow

1. **QC** (optional): CheckM2 subprocess runner for completeness/contamination scoring
1.5. **Re-annotate** (optional): Bakta re-annotation of raw genome assemblies
2. **Parse**: Memory-mapped GFF3/FASTA parsing with Rayon parallelism
3. **Cluster**: MMseqs2-GPU or CPU SIMD fallback
4. **Build Graph**: DashMap concurrent graph construction
5. **Correct**: Paralog resolution, contamination removal, fragment merging, missing gene recovery
6. **Output**: Parallel output generation (CSV, FASTA, GML, JSON, Parquet, HTML)

## Key Data Types

- `Gene`: Parsed gene with sequence, coordinates, and metadata
- `GeneCluster`: Cluster of orthologous genes with centroid sequence
- `PangenomeGraph`: Nodes (DashMap) + edges for the pangenome
- `ConcurrentGraph`: Lock-free concurrent graph using DashMap
- `BitPackedMatrix`: Compressed presence/absence matrix
- `BaktaRunner`: Subprocess runner for Bakta genome annotation
- `BaktaDbType`: Database type selection (Full/Light)

## Configuration Pattern

Configuration uses a builder pattern via `PanminerConfig::new().with_*(...)` methods. See `src/config.rs` for all options.

## Error Handling

Errors use `thiserror` with `Error` enum and `Result<T>` alias. Some errors are "recoverable" (check with `error.is_recoverable()`).

## Testing

PanMiner uses a combination of unit tests and integration tests.

### Unit Tests

Unit tests are integrated directly in the source files under `#[cfg(test)]` blocks. Use `cargo test` to run all unit tests.

### Integration Tests

Integration tests are located in the `tests/` directory. Currently contains a placeholder test — needs expansion with real pipeline test cases.

### Test Status

All tests pass: **117 unit tests + 5 fragment merge + 9 integration + 4 P1 integration + 5 doc tests**

## Priority Roadmap

See [Comparison.md](Comparison.md#priority-roadmap) for the current development priorities:
- **P0 (Critical)**: All completed — fragment merger with real sequences, missing gene recovery with semi-global alignment, full contig DNA, paralog resolution
- **P1 (High)**: Contig-end pruning ✅, structural variant matrix ✅, Mash wrapper for MDS projections, integration tests ✅
- **P2 (Medium)**: GWAS integration, alignment tool integration, Parquet output, HTML visualization
- **P3 (Nice-to-have)**: SIMD optimization, GPU compute, Python bindings, evolutionary models

**Completed in v0.1.0**: CheckM2 pre-processing QC integration, MSA integration (MAFFT/Clustal/PRANK), structural variant detection

**Completed in v0.2.0**: Levenshtein alignment in FragmentMerger, semi-global alignment in MissingGeneRecoverer, full contig DNA from GFF FASTA, paralog resolution with synteny, BFS depth [1,2,3] in mistranslation correction, DistanceCache for reuse, pre-clustering by length/prefix, misassembly edge cleaning, summary statistics

**Completed in v0.3.0**: Bakta re-annotation integration (Phase 0.5), GenBank-to-FASTA conversion, mixed input handling (GFF/FASTA/GenBank)

## Panaroo Reference

PanMiner is inspired by [Panaroo](https://github.com/gtonkinhill/panaroo) (Tonkin-Hill et al., Genome Biology 2020) but uses modern alternatives:

| Aspect | Panaroo | PanMiner |
|--------|---------|----------|
| Language | Python | Rust |
| Clustering | CD-HIT | MMseqs2 + CPU fallback |
| Graph | NetworkX | DashMap + petgraph |
| I/O | Standard reads | Memory-mapped (mmap2) |
| GPU | None | MMseqs2-GPU support |
| Parallelism | multiprocessing | Rayon work-stealing |

See [PANAROO_ALGORITHM.md](PANAROO_ALGORITHM.md) for detailed Panaroo algorithm documentation.