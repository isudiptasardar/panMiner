# PanMiner

A modern pangenome analysis tool written in Rust. PanMiner processes GFF3-annotated genome assemblies to build pangenome graphs, with support for GPU-accelerated clustering via MMseqs2 and CPU fallback.

[![crates.io](https://img.shields.io/crates/v/panminer)](https://crates.io/crates/panminer)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
[![Rust 1.70+](https://img.shields.io/badge/rust-1.70+-lightgrey)](https://www.rust-lang.org)

## Features

- **Memory-mapped I/O** - Zero-copy file access for large datasets
- **Parallel processing** - Rayon-based work-stealing parallelism
- **GPU-accelerated clustering** - MMseqs2 with CUDA support
- **Concurrent graph** - Lock-free DashMap-based graph construction
- **Multiple output formats** - CSV, FASTA, GML, JSON, JSONL, Parquet, HTML

## Installation

### Prerequisites

- **Rust 1.70+** - [Install via rustup](https://rustup.rs)
- **MMseqs2** (optional) - For GPU-accelerated clustering

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
git clone https://github.com/panminer/panminer.git
cd panminer

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

**MMseqs2** (for GPU-accelerated clustering):

```bash
# Via conda
conda install -c bioconda mmseqs2

# Or download from https://github.com/soedinglab/MMseqs2
```

**Note**: PanMiner automatically detects MMseqs2 and NVIDIA GPUs. If MMseqs2 is not installed, it falls back to CPU-based greedy clustering. Use `--no-gpu` or `--no-mmseqs2` to disable GPU/MMseqs2 features.

## Quick Start

```bash
# Basic usage with default settings
panminer genome1.gff genome2.gff -o output_dir

# With custom identity threshold and threads
panminer *.gff -o panminer_output --identity 0.95 --threads 8

# Force CPU mode (disable GPU)
panminer *.gff -o output --force-cpu
```

### Command Line Options

| Flag | Description | Default |
|------|-------------|---------|
| `INPUT` | Input GFF3 files (required) | - |
| `-o, --output` | Output directory | `panminer_output` |
| `-t, --threads` | Thread count (0 = auto) | `0` |
| `--chunk-size` | Genomes per chunk for streaming | `100` |
| `--identity` | Clustering identity threshold (0.5-1.0) | `0.98` |
| `--mode` | Correction mode: strict, default, sensitive | `default` |
| `--force-cpu` | Disable GPU/MMseqs2 | `false` |
| `--no-mmseqs2` | Disable MMseqs2 clustering | `false` |
| `--no-gpu` | Disable GPU detection and acceleration | `false` |
| `--mmseqs-path` | Path to MMseqs2 binary | Auto-detect |
| `-v, --verbose` | Enable debug logging | `false` |
| `--formats` | Output formats (comma-separated) | `matrix,alignment,graph` |

### Correction Modes

| Mode | Contamination Threshold | Use Case |
|------|------------------------|----------|
| `strict` | 5% of genome count | Phylogenetic studies |
| `default` | 2 genomes | Most use cases |
| `sensitive` | 1 genome | Preserve rare genes |

## Output Files

| File | Description |
|------|-------------|
| `_gene_presence_absence.csv` | Gene presence/absence matrix (Roary-compatible) |
| `_core_alignment.fasta` | Core genome alignment |
| `_pangenome_graph.gml` | Cytoscape-compatible graph |
| `_pangenome.json` | Full pangenome summary |
| `_pangenome.jsonl` | Streaming JSON output |

## Usage Examples

### Process multiple genomes with custom settings

```bash
panminer \
  genome1.gff genome2.gff genome3.gff \
  -o panminer_results \
  --identity 0.95 \
  --mode strict \
  --threads 16
```

### Use GPU acceleration (requires MMseqs2 with CUDA)

```bash
panminer *.gff -o gpu_output --mode sensitive
# MMseqs2 GPU is auto-detected and used if available
```

### Stream large datasets (exceeds memory)

```bash
panminer *.gff -o streaming_output --chunk-size 50
```

### Generate specific output formats

```bash
panminer *.gff -o output --formats matrix,json
```

## Library Usage

```rust
use panminer::{PanminerConfig, PanminerPipeline};

fn main() -> panminer::Result<()> {
    let config = PanminerConfig::new()
        .with_input_files(vec!["genome1.gff".into(), "genome2.gff".into()])
        .with_output_dir("output".into())
        .with_threads(8)
        .with_identity(0.98);
    
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
├── io/           # I/O & Memory (mmap, gff, fasta, compress, streaming)
├── clustering/   # Compute - Clustering (mmseqs, cpu, traits)
├── graph/        # Data Structures (types, concurrent, matrix, builder)
├── correction/   # Error Correction (contamination, fragment, missing)
├── output/       # Output Generation (matrix, alignment, graph, json)
└── pipeline.rs   # Main pipeline orchestration
```

See [Architecture.md](Architecture.md) for detailed module descriptions.

## Configuration Reference

See [Specs.md](Specs.md) for:
- Complete feature matrix
- Configuration options
- Correction mode details
- Output format specifications

## Feature Comparison

See [Comparison.md](Comparison.md) for detailed comparison with Panaroo and the PanMiner roadmap.

## Known Gaps

See [Comparison.md](Comparison.md#priority-roadmap) for current limitations:
- Pre-processing QC (Mash, CheckM wrappers)
- Real multiple sequence alignment output
- Downstream analysis (GWAS, evolutionary models)
- Complete integration tests

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
