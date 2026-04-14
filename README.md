# PanMiner

A high-performance pangenome analysis tool written in Rust. PanMiner processes genome assemblies (GFF3, FASTA, or GenBank) to build pangenome graphs with GPU-accelerated clustering, pre-processing QC, error correction, and rich downstream analysis — matching and exceeding Panaroo's feature set.

[![crates.io](https://img.shields.io/crates/v/panminer)](https://crates.io/crates/panminer)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
[![Rust 1.70+](https://img.shields.io/badge/rust-1.70+-lightgrey)](https://www.rust-lang.org)
[![tests](https://img.shields.io/badge/tests-326%20passing-brightgreen)]()

## Features

### Core Pipeline
- **GPU-accelerated clustering** — MMseqs2 with CUDA support, built-in CPU fallback
- **Memory-mapped I/O** — Zero-copy file access for large datasets
- **Concurrent graph** — Lock-free DashMap-based graph construction
- **Streaming pipeline** — Chunked bincode+zstd for datasets larger than RAM
- **Mixed input** — Accept GFF3, FASTA (.fna/.fa), and GenBank (.gb/.gbk/.gbff) files

### Error Correction (6 modules, all wired)
- **Paralog resolution** — BFS context vector similarity (depth 5) + centroid length check
- **Contamination removal** — Iterative low-support degree-1 node removal
- **Contig-end pruning** — Removes terminal low-support boundary nodes
- **Fragment merging** — Mistranslation correction + gene family collapsing with DistanceCache
- **Missing gene recovery** — Semi-global HW alignment in flanking contig sequences
- **Misassembly edge cleaning** — Two-criteria removal (contig-end + disproportionate edges)
- **Highly variable gene detection** — Cycle-based graph analysis (Panaroo-compatible algorithm)

### Pre-processing & QC
- **CheckM2 integration** — Assembly completeness and contamination scoring
- **Distance estimation** — skani (sparse k-mer chaining, 50× faster than FastANI)
- **MDS projection** — Pure Rust classical MDS for genome distance visualization
- **QC visualization** — HTML report with d3.js MDS scatter + bar charts
- **Bakta re-annotation** — Annotate raw genome assemblies before analysis

### Output Formats
- **Roary-compatible CSV** — `gene_presence_absence.csv` (14 metadata columns)
- **Roary gene member CSV** — `gene_presence_absence_roary.csv` (semicolon gene IDs)
- **Binary matrix** — `gene_presence_absence.Rtab`
- **Enriched GML** — Graph with length, seq, protein, genome_ids, member, is_paralog, is_highly_variable
- **Panaroo reference files** — `pan_genome_reference.fa`, `gene_data.csv` (with DNA/protein + location)
- **Combined FASTA** — `combined_DNA_CDS.fasta`, `combined_protein_CDS.fasta`
- **Core genome alignment** — MAFFT, Clustal Omega, or PRANK
- **Alignment trimming** — ClipKIT + BMGE filtering
- **Codon alignment** — MACSE v2 via Java subprocess
- **Structural variant matrix** — `struct_presence_absence.csv` and `.tsv`
- **JSON/JSONL** — `_pangenome.json` and `_pangenome.jsonl`
- **Parquet** — Apache Arrow/Parquet (feature-gated)
- **HTML visualization** — d3.js force-directed graph (feature-gated)
- **Summary statistics** — Core/Soft core/Shell/Cloud classification + highly variable gene count

### Downstream Analysis (`panminer analyze`)
- **Pyseer** — Pan-GWAS via subprocess
- **Scoary2** — Gene-trait association (Genome Biology 2024)
- **SpydrPick** — MI-based epistasis detection (NAR 2019)
- **Panstripe** — Phylogeny-aware gene gain/loss rates (Genome Biology 2023)
- **AMRFinderPlus** — Curated AMR detection (NCBI)
- **Gene neighborhood** — Native BFS extraction from pangenome graph
- **Accumulation curves** — Native rarefaction with Heaps' law fitting
- **GrapeTree/iTOL** — Native profile/annotation export

### Infrastructure
- **Subprocess timeout** — All external tools run with 1-hour timeout protection
- **Error on feature absence** — cDBG mode returns `Error::FeatureNotEnabled` instead of silently continuing
- **cDBG pipeline mode** — GGCAT + ggCaller for de novo gene calling

## Supported Platforms

**Linux** and **macOS** are the primary supported platforms. This aligns with the bioinformatics ecosystem — MMseqs2, CD-HIT, Bakta, CheckM2, skani, and other external tools are available via conda on Linux/macOS only.

Windows users should use **WSL2** (Windows Subsystem for Linux) or run PanMiner on a Linux server/HPC cluster.

## Installation

### Prerequisites

- **Rust 1.70+** — [Install via rustup](https://rustup.rs)
- **MMseqs2** (optional) — For GPU-accelerated clustering (`conda install -c bioconda mmseqs2`)
- **skani** (optional) — For fast, robust ANI distance estimation (`conda install -c bioconda skani`)
- **Bakta** (optional) — For re-annotation of raw genome assemblies
- **CheckM2** (optional) — For pre-processing quality control
- **MAFFT/Clustal/PRANK** (optional) — For multiple sequence alignment
- **ClipKIT** (optional) — For alignment trimming
- **pyseer** (optional) — For pan-GWAS analysis

### Method 1: Install from Cargo

```bash
cargo install panminer
```

### Method 2: Clone and Build

```bash
git clone https://github.com/isudiptasardar/panMiner.git
cd panMiner
cargo build --release
cargo install --path .
```

### Method 3: Conda Environment (Recommended for Full Toolchain)

**Rust must be installed first** (via rustup, not conda) before creating the conda environment:

```bash
# Step 1: Install Rust (required before conda env)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"

# Step 2: Create conda environment with external tools
conda env create -f environment.yml
conda activate panminer

# Step 3: Build and install PanMiner
cargo build --release
cargo install --path .
```

The `environment.yml` installs PanMiner's external tool dependencies (MMseqs2, skani, MAFFT, etc.) but **not Rust itself** — Rust should come from rustup to avoid conda solver conflicts.

**Troubleshooting:** If conda takes too long or hangs during solving:
1. Install with `conda env create -f environment.yml --strict-channel-priority`
2. If that fails, comment out more packages in `environment.yml` and retry
3. Alternatively, install individual tools as needed (see External Tools below)

**Risky packages** (checkm2, bakta, ggcaller, grapetree) are commented out in `environment.yml` because they often cause conda solver conflicts. Uncomment them individually if needed, or install them separately.

### Method 4: Installation Script (Auto-detect GPU)

```bash
bash install.sh              # Auto-detect GPU, offer MMseqs2 with CUDA
bash install.sh --no-gpu     # Skip GPU detection
bash install.sh --dev        # Build in debug mode
bash install.sh --uninstall  # Remove PanMiner
```

The installation script:
- Checks for Rust and installs if needed
- Detects NVIDIA GPUs automatically
- Prompts to install MMseqs2 with GPU support if a GPU is found
- Builds and installs PanMiner

### External Tools (Optional)

All external tools are auto-detected. If a tool is not installed, its feature is gracefully skipped.

The easiest approach is to use the conda environment (Method 3 above). To install tools individually:

```bash
# Core tools (included in environment.yml)
conda install -c conda-forge -c bioconda mmseqs2 skani mafft clustalo prank clipkit

# Alignment filtering
pip install bmge                                  # BMGE (requires Biopython)
pip install biopython                              # Biopython dependency for BMGE

# Downstream analysis
conda install -c conda-forge -c bioconda pyseer   # Pan-GWAS
conda install -c conda-forge -c bioconda scoary2  # Gene-trait association
conda install -c conda-forge -c bioconda spydrpick # Epistasis detection
conda install -c conda-forge -c bioconda r-base=4.3 r-panstripe  # Evolutionary model
conda install -c conda-forge -c bioconda ncbi-amrfinder  # AMR detection

# Tools with large dep trees — install separately if needed
conda install -c conda-forge -c bioconda checkm2  # QC (large dep tree)
conda install -c conda-forge -c bioconda bakta    # Re-annotation (~6GB DB)
bakta_db download --output ~/.bakta --type full   # Bakta database download
```

## Quick Start

### Basic pangenome analysis

```bash
panminer genome1.gff genome2.gff genome3.gff -o results
```

### Re-annotate raw assemblies

```bash
# Annotate FASTA/GenBank files with Bakta before analysis
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

### With QC and distance estimation

```bash
panminer *.gff -o results --qc-mode default
# Uses skani (if installed) for ANI distance estimation
```

### GPU-accelerated clustering

```bash
# MMseqs2 GPU is auto-detected and used if available
panminer *.gff -o gpu_output --mode sensitive
```

### Stream large datasets

```bash
panminer *.gff -o streaming_output --chunk-size 50
```

### cDBG pipeline mode (de novo gene calling)

```bash
# Requires --features dbg build flag
panminer *.fasta -o dbg_output --pipeline-mode dbg --kmer-size 31
```

### Generate specific output formats

```bash
panminer *.gff -o output --formats matrix,json,graph
```

### Alignment trimming and filtering

```bash
# ClipKIT trimming
panminer *.gff -o output --trim-alignment --trim-mode smart-gap

# BMGE filtering
panminer *.gff -o output --filter-alignment bmge

# Codon alignment (MACSE)
panminer *.gff -o output --codons
```

## Downstream Analysis (`panminer analyze`)

After running `panminer` to build a pangenome, use `panminer analyze` for downstream analysis:

```bash
# Gene neighborhood extraction
panminer analyze -i results --neighborhood --seed-gene cluster_0042 --neighborhood-depth 3

# Accumulation curves
panminer analyze -i results --accumulation --num-samples 100

# GrapeTree/iTOL export
panminer analyze -i results --export-grapetree --export-itol

# GWAS with pyseer
panminer analyze -i results --gwas --gwas-tool pyseer --phenotypes phenotypes.tsv

# Gene-trait association with Scoary2
panminer analyze -i results --gwas --gwas-tool scoary2 --phenotypes phenotypes.tsv

# Epistasis detection with SpydrPick
panminer analyze -i results --gwas --gwas-tool spydrpick --phenotypes phenotypes.tsv

# Evolutionary model with Panstripe (requires phylogenetic tree)
panminer analyze -i results --panstripe --tree tree.nwk

# AMR detection with AMRFinderPlus
panminer analyze -i results --amr --organism "Escherichia coli"

# P/A filtering
panminer filter-pa -i results -o filtered_output
```

## Command Line Reference

### Main Pipeline (`panminer`)

| Flag | Description | Default |
|------|-------------|---------|
| `INPUT` | Input GFF3/FASTA/GenBank files (required) | — |
| `-o, --output` | Output directory | `panminer_output` |
| `-t, --threads` | Thread count (0 = auto) | `0` |
| `--chunk-size` | Genomes per chunk for streaming | `100` |
| `--identity` | Clustering identity threshold (0.5–1.0) | `0.98` |
| `--mode` | Correction mode: strict, default, sensitive | `default` |
| `--pipeline-mode` | Pipeline mode: gff or dbg | `gff` |
| `--kmer-size` | K-mer size for cDBG mode (15–127) | `31` |
| `-r, --reannotate` | Re-annotate inputs with Bakta | `false` |
| `--bakta-db` | Path to Bakta database directory | Auto-detect |
| `--bakta-db-type` | Bakta DB type: full or light | `full` |
| `--bakta-threads` | Threads for Bakta | Auto |
| `--no-bakta-db-download` | Fail if Bakta DB not found | `false` |
| `--keep-bakta-output` | Keep Bakta temp files | `false` |
| `--force-cpu` | Disable GPU/MMseqs2 | `false` |
| `--no-mmseqs2` | Disable MMseqs2 clustering | `false` |
| `--no-gpu` | Disable GPU detection | `false` |
| `--no-qc` | Disable pre-processing QC | `false` |
| `--qc-mode` | QC mode: strict, default, sensitive | `default` |
| `--checkm-database` | Path to CheckM2 database | Auto-detect |
| `--mmseqs-path` | Path to MMseqs2 binary | Auto-detect |
| `--verbose` | Enable debug logging | `false` |
| `--formats` | Output formats (comma-separated) | `matrix,alignment,graph` |
| `--trim-alignment` | Trim alignment with ClipKIT | `false` |
| `--trim-mode` | ClipKIT trim mode | `smart-gap` |
| `--filter-alignment` | Alignment filter: none, clipkit, bmge | `none` |
| `--codons` | Generate codon alignments (MACSE) | `false` |
| `--gwas` | Run GWAS with pyseer | `false` |
| `--phenotype` | Phenotype file for GWAS | None |

### Downstream Analysis (`panminer analyze`)

| Flag | Description |
|------|-------------|
| `-i, --input <dir>` | PanMiner output directory (required) |
| `--neighborhood` | Extract gene neighborhood |
| `--seed-gene <id>` | Seed cluster ID for neighborhood extraction |
| `--neighborhood-depth <n>` | Maximum BFS depth (default: 5) |
| `--accumulation` | Generate gene accumulation curves |
| `--num-samples <n>` | Rarefaction samples per point (default: 100) |
| `--export-grapetree` | Export GrapeTree profiles |
| `--export-itol` | Export iTOL annotations |
| `--gwas` | Run GWAS analysis |
| `--gwas-tool <tool>` | GWAS tool: pyseer, scoary2, spydrpick |
| `--phenotypes <file>` | Phenotypes TSV (genome_id tab phenotype) |
| `--panstripe` | Run Panstripe evolutionary model |
| `--tree <file>` | Phylogenetic tree (Newick format) |
| `--amr` | Run AMRFinderPlus resistome analysis |
| `--amr-database <path>` | AMRFinderPlus database path |
| `--organism <name>` | Organism for taxon-specific AMR |

### P/A Filtering (`panminer filter-pa`)

| Flag | Description |
|------|-------------|
| `-i, --input <dir>` | PanMiner output directory (required) |
| `-o, --output <dir>` | Filtered output directory |
| `--min-genomes <n>` | Minimum number of genomes a gene must be in |
| `--max-genomes <n>` | Maximum number of genomes a gene must be in |
| `--remove-fragments` | Remove fragmented genes |
| `--remove-pseudogenes` | Remove pseudogenes |

## Correction Modes

| Mode | Contamination Threshold | Use Case |
|------|------------------------|----------|
| `strict` | 5% of genome count | Phylogenetic studies — removes more |
| `default` | 2 genomes | Most use cases — balanced |
| `sensitive` | 1 genome | Preserve rare genes/plasmids |

## Distance Estimation Priority

PanMiner uses **skani** for distance estimation — the fastest and most robust ANI tool available.

1. **skani** — 50× faster than FastANI, robust to incomplete genomes (MAGs)

Install: `conda install -c bioconda skani`

If skani is not installed, the QC distance step is skipped.

## Output Files

PanMiner follows **Panaroo/Roary naming conventions** for compatibility with downstream tools:

### Core Output Files

| File | Description | Roary Compatible |
|------|-------------|------------------|
| `gene_presence_absence.csv` | Gene P/A matrix (14 metadata columns) | Yes |
| `gene_presence_absence_roary.csv` | Gene P/A with semicolon gene IDs per genome | Yes |
| `gene_presence_absence.Rtab` | Binary tab-separated presence/absence | Yes |
| `final_graph.gml` | Pangenome graph (enriched GML) | Yes |
| `pre_filt_graph.gml` | Graph before correction | Yes |
| `core_gene_alignment.aln` | Core gene alignment (MSA) | Yes |
| `core_gene_alignment.trimmed.aln` | ClipKIT-trimmed alignment | No |
| `core_gene_alignment.BMGE.aln` | BMGE-filtered alignment | No |
| `core_gene_alignment.codon.aln` | MACSE codon alignment | No |
| `pan_genome_reference.fa` | Reference genome of all genes | Yes |
| `gene_data.csv` | Gene-to-cluster mapping + DNA/protein sequences | No |
| `combined_DNA_CDS.fasta` | All nucleotide sequences | No |
| `combined_protein_CDS.fasta` | All protein sequences | No |
| `struct_presence_absence.csv` | Gene triplet structural variants | No |
| `struct_presence_absence.Rtab` | Binary structural variant matrix | No |
| `summary_statistics.txt` | Core/Soft core/Shell/Cloud + highly variable counts | No |

### QC Output Files (when enabled)

| File | Description |
|------|-------------|
| `qc_stats.csv` | Per-genome QC metrics (completeness, contamination, ANI, etc.) |
| `qc_summary.txt` | Human-readable QC summary with pass/fail status |
| `qc_viz.html` | Interactive HTML report with MDS scatter + bar charts |

### Additional Files

| File | Description |
|------|-------------|
| `_pangenome.json` | Full pangenome summary (JSON) |
| `_pangenome.jsonl` | Streaming JSON output (one cluster per line) |
| `_pangenome.parquet` | Parquet format (requires `--features parquet`) |
| `_pangenome.html` | Interactive HTML visualization (requires `--features viz`) |

### GML Graph Attributes

The `final_graph.gml` file includes these node attributes:

| Attribute | Description |
|-----------|-------------|
| `id` / `label` | Cluster ID |
| `support` | Number of genomes containing this gene |
| `is_paralog` | Whether this is a paralogous cluster (0/1) |
| `is_highly_variable` | Whether this is a highly variable gene (0/1) |
| `length` | Length of centroid sequence |
| `seq` | Centroid DNA sequence |
| `protein` | Centroid protein sequence |
| `genome_ids` | Genomes containing this gene (GML list) |
| `member` | Per-genome gene IDs (GML list) |

## Pipeline Flow

```
Phase 0:    QC (optional) — CheckM2 completeness/contamination + skani ANI distance
Phase 0.5:  Re-annotation (optional) — Bakta annotation of raw assemblies
Phase 1:    Parse — Memory-mapped GFF3/FASTA parsing with Rayon parallelism
Phase 2:    Cluster — MMseqs2-GPU or CPU fallback greedy clustering
Phase 3:    Build graph — DashMap concurrent graph construction
Phase 4:    Correct — 6 correction modules:
             4.0 Paralog resolution (BFS context vectors)
             4.1 Contamination removal
             4.2 Contig-end pruning
             4.3 Fragment merging (mistranslation + gene family collapse)
             4.4 Missing gene recovery (semi-global HW alignment)
             4.5 Re-collapse gene families
             4.6 Misassembly edge cleaning
Phase 4.8:  Detect highly variable genes (cycle-based graph analysis)
Phase 5:    Matrix — Build presence/absence matrix (BitPackedMatrix)
Phase 6:    Output — Generate all output files
Phase 7:    GWAS (optional) — Pyseer/Scoary2/SpydrPick subprocess
```

## Architecture

```
src/
├── lib.rs              # Library entry point, re-exports
├── main.rs             # CLI with clap
├── config.rs           # Configuration structs
├── error.rs            # Error types (thiserror-based)
├── pipeline.rs         # Main pipeline orchestration
│
├── io/                 # I/O & Memory
│   ├── mmap.rs         # Memory-mapped file wrapper
│   ├── gff.rs          # GFF3 parser (mmap-based)
│   ├── fasta.rs        # FASTA parser
│   ├── compress.rs     # Zstd compression
│   ├── streaming.rs    # Chunked processing for large datasets
│   ├── subprocess.rs   # Centralized subprocess execution with timeout
│   ├── bakta.rs        # Bakta re-annotation runner
│   ├── qc_traits.rs    # QC runner traits, GenomeQC, ANI dispatch
│   ├── translate.rs    # DNA-to-protein translation
│   ├── skani.rs        # skani ANI subprocess (sole distance tool)
│   ├── mds.rs          # Classical MDS projection (pure Rust)
│   ├── ggcat.rs        # GGCAT colored cDBG (feature-gated)
│   └── ggcaller.rs     # ggCaller gene calling subprocess
│
├── clustering/         # Compute — Clustering
│   ├── traits.rs       # Clusterer trait
│   ├── mmseqs.rs       # MMseqs2-GPU integration
│   ├── cpu.rs          # CPU fallback with SIMD
│   └── alignment_*.rs  # MSA tools (MAFFT, Clustal, PRANK)
│
├── graph/              # Data Structures
│   ├── types.rs        # Gene, GeneCluster, Node, Edge, PangenomeGraph
│   ├── concurrent.rs   # DashMap-based ConcurrentGraph
│   ├── matrix.rs       # BitPackedMatrix (8x memory reduction)
│   ├── builder.rs      # Graph construction with gene_members
│   ├── structural_variants.rs  # Structural variant detection
│   ├── highly_variable.rs     # Highly variable gene detection
│   └── merge.rs        # Pangenome merging
│
├── correction/         # Error Correction (6 modules)
│   ├── contamination.rs    # Low-support node removal
│   ├── contig_end.rs       # Contig-end pruning
│   ├── fragment.rs        # Fragment merging with DistanceCache
│   ├── missing.rs          # Missing gene recovery (semi-global HW)
│   ├── misassembly.rs     # Misassembly edge cleaning
│   ├── paralog.rs          # Paralog resolution (BFS context vectors)
│   └── simd.rs             # SIMD sequence comparison + Levenshtein
│
├── downstream/         # Downstream Analysis
│   ├── gwas/            # PyseerRunner, Scoary2Runner, SpydrPickRunner
│   ├── evolution/       # PanstripeRunner
│   ├── resistome/       # AmrFinderRunner
│   └── exploration/     # AccumulationCurve, GeneNeighborhood, GrapeTree, iTOL
│
└── output/              # Output Generation
    ├── matrix.rs        # Roary CSV (14-col), Rtab, gene member CSV
    ├── alignment.rs     # Core/accessory alignment via MSA
    ├── graph.rs         # GML (enriched attributes)
    ├── json.rs          # gene_data.csv, pan_genome_reference.fa, JSON
    ├── struct_csv.rs    # Structural variant CSV
    ├── sv_matrix.rs     # Structural variant TSV
    ├── summary.rs       # Core/Soft core/Shell/Cloud + highly variable
    ├── parquet.rs       # Apache Parquet (feature-gated)
    ├── html_viz.rs      # d3.js HTML visualization (feature-gated)
    ├── filter_pa.rs     # P/A filtering
    ├── trim.rs          # ClipKIT + BMGE filtering
    ├── codon.rs         # MACSE codon alignment
    ├── qc_stats.rs      # QC statistics
    └── qc_viz.rs        # QC HTML report
```

## Panaroo Comparison

| Aspect | Panaroo | PanMiner |
|--------|---------|----------|
| Language | Python | Rust |
| Clustering | CD-HIT | MMseqs2 (GPU) + CPU fallback |
| Graph | NetworkX | DashMap concurrent + petgraph |
| Distance | Mash | skani (sparse k-mer chaining) |
| I/O | Standard reads | Memory-mapped (mmap2) |
| Memory | High (NetworkX) | Low (BitPackedMatrix 8× reduction) |
| Streaming | No | Yes (chunked bincode+zstd) |
| GPU | None | MMseqs2-GPU subprocess |
| Parallelism | multiprocessing | Rayon work-stealing |
| Pipeline modes | GFF3 only | GFF3 + cDBG (GGCAT + ggCaller) |
| Error correction | 6 modules | 6 modules + highly variable detection |
| Downstream tools | Scoary, SpydrPick | Scoary2, SpydrPick, Panstripe, AMRFinderPlus |

## Build Features

| Feature | Description |
|---------|-------------|
| `cpu` (default) | CPU processing support |
| `mmseqs` | MMseqs2 integration for clustering |
| `parquet` | Parquet output format support |
| `python` | PyO3 Python bindings |
| `viz` | HTML visualization output |
| `dbg` | GGCAT cDBG pipeline mode |
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
cargo test test_highly_variable

# Run tests with all features
cargo test --features full

# Check without building
cargo check

# Run lints
cargo clippy

# Format code
cargo fmt
```

### Documentation

```bash
# Build and open documentation
cargo doc --open
```

See [CLAUDE.md](CLAUDE.md), [Specs.md](Specs.md), [Architecture.md](Architecture.md), and [Comparison.md](Comparison.md) for detailed development documentation.

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

MIT License — see [LICENSE](LICENSE) for details.

## Acknowledgments

- PanMiner is inspired by [Panaroo](https://github.com/gtonkinhill/panaroo) (Tonkin-Hill et al., Genome Biology 2020)
- Built with Rust, Rayon, DashMap, memmap2, and petgraph
- Uses MMseqs2 for GPU-accelerated clustering (Steinegger & Söding, Nature Biotechnology 2017)
- Uses skani for fast ANI estimation (Shaw & Yu, Nature Methods 2023)
- Uses Bakta for genome re-annotation (Schwengers et al., Microbial Genomics 2021)