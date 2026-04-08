# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PanMiner (panminer) is a modern pangenome analysis tool written in Rust. It processes GFF3 files to build pangenome graphs, with support for GPU-accelerated clustering via MMseqs2 and CPU fallback. The project follows a layered architecture inspired by Panaroo but with modern optimizations (memory-mapped I/O, DashMap concurrent graphs, Rayon parallelism).

## Build and Test Commands

```bash
# Build
cargo build

# Build release (optimized)
cargo build --release

# Run tests
cargo test

# Run specific test
cargo test test_name

# Run tests with features
cargo test --features full

# Check without building
cargo check

# Run the binary with help
cargo run -- --help

# Run with specific inputs
cargo run -- input1.gff input2.gff -o output_dir
```

## Feature Flags

- `cpu` (default): CPU processing support
- `mmseqs`: MMseqs2 integration for clustering
- `parquet`: Parquet output format support
- `python`: PyO3 Python bindings
- `viz`: HTML visualization output
- `full`: All features enabled

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
│   └── streaming.rs    # Chunked processing for large datasets
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
│   └── missing.rs
│
└── output/             # Layer 4: Output Generation
    ├── matrix.rs       # CSV/TSV presence-absence
    ├── alignment.rs    # FASTA alignments
    ├── graph.rs        # GML output
    └── json.rs         # JSON/JSONL output
```

## Pipeline Flow

1. **Parse**: Memory-mapped GFF3/FASTA parsing with Rayon parallelism
2. **Cluster**: MMseqs2-GPU or CPU SIMD fallback
3. **Build Graph**: DashMap concurrent graph construction
4. **Correct**: Contamination removal, fragment merging, missing gene recovery
5. **Output**: Parallel output generation (CSV, FASTA, GML, JSON, Parquet, HTML)

## Key Data Types

- `Gene`: Parsed gene with sequence, coordinates, and metadata
- `GeneCluster`: Cluster of orthologous genes with centroid sequence
- `PangenomeGraph`: Nodes (DashMap) + edges for the pangenome
- `ConcurrentGraph`: Lock-free concurrent graph using DashMap
- `BitPackedMatrix`: Compressed presence/absence matrix

## Configuration Pattern

Configuration uses a builder pattern via `PanminerConfig::new().with_*(...)` methods. See `src/config.rs` for all options.

## Error Handling

Errors use `thiserror` with `Error` enum and `Result<T>` alias. Some errors are "recoverable" (check with `error.is_recoverable()`).

## Testing

Tests are inline in module files under `#[cfg(test)]` blocks. Use `cargo test` to run all tests. Development test scripts have been moved to the `scripts/` directory (`scripts/add_*_test.sh`, `scripts/fix_*.sh`) and can be used during development to add tests incrementally.