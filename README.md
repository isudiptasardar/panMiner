# PanMiner

A modern pangenome analysis tool written in Rust. PanMiner processes genome assemblies (GFF3, FASTA, or GenBank) to build pangenome graphs, with support for GPU-accelerated clustering via MMseqs2, pre-processing QC via CheckM2, re-annotation via Bakta, and CPU fallback.

[![crates.io](https://img.shields.io/crates/v/panminer)](https://crates.io/crates/panminer)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
[![Rust 1.70+](https://img.shields.io/badge/rust-1.70+-lightgrey)](https://www.rust-lang.org)

## Features

- **Bakta re-annotation** - Annotate raw genome assemblies (FASTA/GenBank) before analysis
- **Pre-processing QC** - CheckM2 integration for completeness/contamination scoring
- **Memory-mapped I/O** - Zero-copy file access for large datasets
- **Parallel processing** - Rayon-based work-stealing parallelism
- **GPU-accelerated clustering** - MMseqs2 with CUDA support, built-in CPU fallback
- **Concurrent graph** - Lock-free DashMap-based graph construction
- **Error correction** - Paralog resolution, contig-end pruning, fragment merging, missing gene recovery, misassembly edge cleaning
- **Structural variant matrix** - Gene triplet presence/absence output
- **Multiple output formats** - CSV, FASTA, GML, JSON, JSONL, Parquet, HTML
- **Multiple alignment tools** - MAFFT, Clustal Omega, PRANK subprocess integration
- **Mixed input** - Accept GFF3, FASTA (.fna/.fa), and GenBank (.gb/.gbk/.gbff) files

## Installation

### Prerequisites

- **Rust 1.70+** - [Install via rustup](https://rustup.rs)
- **MMseqs2** (optional) - For GPU-accelerated clustering
- **Bakta** (optional) - For re-annotation of raw genome assemblies
- **CheckM2** (optional) - For pre-processing quality control

### Method 1: Install from Cargo (Recommended)

```bash
# Build and install from crates.io
cargo install panminer

# Or build and install from source
cargo install --path .
```

### Method 2: Clone and Build

```bash
# Clone the repository
git clone https://github.com/isudiptasardar/panMiner.git
cd panMiner

# Build in release mode
cargo build --release

# Install to your PATH
cargo install --path .
```

### Method 3: Conda Environment

```bash
# Create conda environment with all dependencies
conda env create -f environment.yml
conda activate panminer

# Build and install
cargo build --release
cargo install --path .
```

### Installation Script (Recommended)

PanMiner includes an installation script that automatically detects your system and installs dependencies:

```bash
# Basic installation
bash install.sh

# With options
bash install.sh --dev           # Build in debug mode
bash install.sh --no-gpu        # Skip GPU detection and MMseqs2 installation
bash install.sh --no-mmseqs2    # Skip MMseqs2 installation entirely
bash install.sh --uninstall     # Remove PanMiner
```

The installation script will:
- Check for Rust and install if needed
- Detect NVIDIA GPUs automatically
- Prompt to install MMseqs2 with GPU support if a GPU is found
- Build and install PanMiner

### External Tool Installation (Optional)

#### CheckM2 (Pre-processing QC)
```bash
# Via conda
conda install -c bioconda checkm2

# Or download from https://github.com/chklovski/CheckM2
```

#### MMseqs2 (GPU-accelerated clustering)
```bash
# Via conda
conda install -c bioconda mmseqs2

# Or download from https://github.com/soedinglab/MMseqs2
```

#### Bakta (Genome re-annotation)
```bash
# Via conda (recommended)
conda install -c conda-forge -c bioconda bakta

# Or via pip
pip install bakta

# Download the Bakta database (full or light)
bakta_db download --output ~/.bakta --type full
```

**Note**: PanMiner automatically detects CheckM2, MMseqs2, and Bakta. If a tool is not installed, its feature is gracefully skipped. Use `--no-qc`, `--no-mmseqs2`, or simply don't pass `--reannotate` to disable specific features.

## Quick Start

```bash
# Basic usage with GFF3 files
panminer genome1.gff genome2.gff -o output_dir

# With custom identity threshold and threads
panminer *.gff -o panminer_output --identity 0.95 --threads 8

# Re-annotate raw assemblies with Bakta before analysis
panminer -r genome1.fasta genome2.fasta genome3.gff -o output_dir

# Force CPU mode (disable GPU)
panminer *.gff -o output --force-cpu
```

### Command Line Options

| Flag | Description | Default |
|------|-------------|---------|
| `INPUT` | Input GFF3/FASTA/GenBank files (required) | - |
| `-o, --output` | Output directory | `panminer_output` |
| `-t, --threads` | Thread count (0 = auto) | `0` |
| `--chunk-size` | Genomes per chunk for streaming | `100` |
| `--identity` | Clustering identity threshold (0.5-1.0) | `0.98` |
| `--mode` | Correction mode: strict, default, sensitive | `default` |
| `-r, --reannotate` | Re-annotate inputs with Bakta before analysis | `false` |
| `--bakta-db` | Path to Bakta database directory | Auto-detect |
| `--bakta-db-type` | Bakta DB type for auto-download: full or light | `full` |
| `--bakta-threads` | Threads for Bakta (default: same as pipeline) | Auto |
| `--no-bakta-db-download` | Fail if Bakta DB not found (no auto-download) | `false` |
| `--keep-bakta-output` | Keep Bakta output files after pipeline | `false` |
| `--force-cpu` | Disable GPU/MMseqs2 | `false` |
| `--no-mmseqs2` | Disable MMseqs2 clustering | `false` |
| `--no-gpu` | Disable GPU detection and acceleration | `false` |
| `--no-qc` | Disable pre-processing QC | `false` |
| `--qc-mode` | QC mode: strict, default, sensitive | `default` |
| `--checkm-database` | Path to CheckM2 database | Auto-detect |
| `--mmseqs-path` | Path to MMseqs2 binary | Auto-detect |
| `-v, --verbose` | Enable debug logging | `false` |
| `--formats` | Output formats (comma-separated) | `matrix,alignment,graph` |

### Re-annotation with Bakta

When `--reannotate` (or `-r`) is passed, PanMiner annotates raw genome assemblies with Bakta before pangenome analysis:

- **GFF/GFF3 files** — passed through unchanged (already annotated)
- **FASTA files** (.fasta, .fna, .fa) — annotated by Bakta
- **GenBank files** (.gb, .gbk, .gbff) — converted to FASTA, then annotated by Bakta

If Bakta is not installed, GFF files are used directly and a warning is logged. GenBank files without Bakta produce an error.

Bakta database resolution priority: `--bakta-db` flag > `BAKTA_DB` env var > `~/.bakta/db` > auto-download (unless `--no-bakta-db-download`).

### Correction Modes

| Mode | Contamination Threshold | Use Case |
|------|------------------------|----------|
| `strict` | 5% of genome count | Phylogenetic studies |
| `default` | 2 genomes | Most use cases |
| `sensitive` | 1 genome | Preserve rare genes |

## Output Files

PanMiner follows **Panaroo/Roary output naming conventions** for compatibility:

| File | Description | Roary Compatible |
|------|-------------|------------------|
| `gene_presence_absence.csv` | Gene presence/absence matrix | Yes |
| `gene_presence_absence.Rtab` | Binary tab-separated presence/absence | Yes |
| `final_graph.gml` | Pan-genome graph (Cytoscape-compatible) | No |
| `pre_filt_graph.gml` | Graph before correction | No |
| `struct_presence_absence.csv` | Genomic rearrangement events | No |
| `struct_presence_absence.Rtab` | Binary structural variant matrix | No |
| `pan_genome_reference.fa` | Reference genome of all genes | Yes (similar) |
| `gene_data.csv` | Gene-to-annotation links | No |
| `combined_DNA_CDS.fasta` | All nucleotide sequences | No |
| `combined_protein_CDS.fasta` | All protein sequences | No |
| `core_gene_alignment.aln` | Core gene alignment | No (uses .aln) |
| `summary_statistics.txt` | Core/Soft core/Shell/Cloud counts | No |

### QC Output Files (when enabled)

| File | Description |
|------|-------------|
| `qc_stats.csv` | Per-genome QC metrics (completeness, contamination, etc.) |
| `qc_summary.txt` | Human-readable QC summary with pass/fail status |

### Additional PanMiner-Specific Files

| File | Description |
|------|-------------|
| `_pangenome.json` | Full pangenome summary (JSON) |
| `_pangenome.jsonl` | Streaming JSON output |
| `_pangenome.parquet` | Parquet format (if `--features parquet`) |
| `_pangenome.html` | Interactive HTML visualization (if `--features viz`) |

## Pipeline Flow

```
Phase 0:   QC (optional) — CheckM2 completeness/contamination scoring
Phase 0.5: Re-annotation (optional) — Bakta annotation of raw assemblies
Phase 1:   Parse — Memory-mapped GFF3/FASTA parsing with Rayon parallelism
Phase 2:   Cluster — MMseqs2-GPU or CPU fallback greedy clustering
Phase 3:   Build graph — DashMap concurrent graph construction
Phase 4:   Correct — Paralog resolution, contamination removal, fragment merging,
                     missing gene recovery, misassembly edge cleaning
Phase 5:   Matrix — Build presence/absence matrix (BitPackedMatrix)
Phase 6:   Output — Generate all output files
```

## Usage Examples

### Basic pangenome analysis

```bash
panminer genome1.gff genome2.gff genome3.gff -o results
```

### Re-annotate raw assemblies with Bakta

```bash
panminer -r *.fasta -o results --threads 16
```

### Mixed input (GFF + FASTA + GenBank)

```bash
panminer -r annotated.gff raw.fasta draft.gbk -o results
```

### Strict mode for phylogenetic studies

```bash
panminer *.gff -o results --mode strict --identity 0.95
```

### GPU acceleration with MMseqs2

```bash
panminer *.gff -o gpu_output --mode sensitive
# MMseqs2 GPU is auto-detected and used if available
```

### Stream large datasets

```bash
panminer *.gff -o streaming_output --chunk-size 50
```

### Generate specific output formats

```bash
panminer *.gff -o output --formats matrix,json,graph
```

## Library Usage

```rust
use panminer::{PanminerConfig, PanminerPipeline, BaktaDbType};
use std::path::PathBuf;

fn main() -> panminer::Result<()> {
    let config = PanminerConfig::new()
        .with_input_files(vec![
            PathBuf::from("genome1.gff"),
            PathBuf::from("genome2.gff"),
        ])
        .with_output_dir(PathBuf::from("output"))
        .with_threads(8)
        .with_reannotate(true)
        .with_bakta_db_type(BaktaDbType::Full);

    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run()?;

    println!("Output: {:?}", result.output_dir);
    Ok(())
}
```

## Build Features

| Feature | Description |
|---------|-------------|
| `cpu` (default) | CPU processing support |
| `mmseqs` | MMseqs2 integration for clustering |
| `parquet` | Parquet output format support |
| `python` | PyO3 Python bindings |
| `viz` | HTML visualization output |
| `full` | All features enabled |

```bash
# Build with specific features
cargo build --features mmseqs,parquet

# Build with all features
cargo build --features full
```

## Development

### Running Tests

```bash
# Run all tests
cargo test

# Run specific test
cargo test test_name

# Run tests with all features
cargo test --features full
```

### Code Quality

```bash
# Check without building
cargo check

# Run lints
cargo clippy

# Format code
cargo fmt
```

### Documentation

```bash
# Build documentation
cargo doc --open
```

## Project Structure

```
src/
├── lib.rs              # Library entry point, re-exports
├── main.rs             # CLI with clap
├── config.rs           # Configuration structs (CorrectionMode, OutputFormat, BaktaDbType)
├── error.rs            # Error types (thiserror-based)
├── pipeline.rs         # Main pipeline orchestration (PanminerPipeline)
│
├── io/                 # I/O & Memory
│   ├── mmap.rs         # Memory-mapped file wrapper
│   ├── gff.rs          # GFF3 parser (mmap-based)
│   ├── fasta.rs        # FASTA parser
│   ├── compress.rs     # Zstd compression utilities
│   ├── streaming.rs    # Chunked processing for large datasets
│   ├── bakta.rs        # Bakta annotation runner (subprocess)
│   └── qc_traits.rs    # QC runner traits (CheckM2)
│
├── clustering/         # Compute
│   ├── traits.rs       # Clusterer trait
│   ├── mmseqs.rs       # MMseqs2-GPU integration
│   ├── cpu.rs          # CPU fallback with SIMD
│   └── alignment_*.rs  # MSA tools (MAFFT, Clustal, PRANK)
│
├── graph/              # Data Structures
│   ├── types.rs        # Gene, GeneCluster, Node, Edge, PangenomeGraph
│   ├── concurrent.rs   # DashMap-based graph (ConcurrentGraph)
│   ├── matrix.rs       # BitPackedMatrix (8x memory reduction)
│   ├── builder.rs      # Graph construction
│   └── structural_variants.rs  # Gene triplet detection
│
├── correction/         # Error Correction
│   ├── contamination.rs    # Low-support node removal
│   ├── contig_end.rs      # Contig-end pruning
│   ├── fragment.rs        # Fragment merging with alignment
│   ├── missing.rs         # Missing gene recovery
│   ├── paralog.rs         # Paralog resolution with synteny
│   ├── misassembly.rs     # Misassembly edge cleaning
│   └── simd.rs            # SIMD sequence comparison
│
└── output/             # Output Generation
    ├── matrix.rs       # CSV/TSV presence-absence
    ├── alignment.rs    # FASTA alignments (MAFFT/Clustal/PRANK)
    ├── graph.rs        # GML output
    ├── json.rs         # JSON/JSONL output
    ├── summary.rs      # Summary statistics
    ├── sv_matrix.rs    # Structural variant matrix
    ├── qc_stats.rs     # QC statistics output
    ├── parquet.rs      # Parquet output
    └── html_viz.rs     # HTML visualization
```

## Known Limitations

- Mash distance estimation (MDS projections) not yet implemented
- GWAS integration (pyseer) not yet implemented
- Python bindings (PyO3) are a stub
- HTML visualization is a stub
- Parquet output is a stub

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- PanMiner is inspired by [Panaroo](https://github.com/gtonkinhill/panaroo)
- Built with Rust, Rayon, DashMap, and memmap2
- Uses MMseqs2 for GPU-accelerated clustering
- Uses Bakta for genome re-annotation