# PanMiner

A high-performance pangenome analysis tool written in Rust. PanMiner processes genome assemblies (GFF3, FASTA, or GenBank) to build pangenome graphs with GPU-accelerated clustering, pre-processing QC, a 7-stage error correction pipeline, and rich downstream analysis -- achieving full Panaroo feature parity with significant architectural improvements.

[![crates.io](https://img.shields.io/crates/v/panminer)](https://crates.io/crates/panminer)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
[![Rust 1.70+](https://img.shields.io/badge/rust-1.70+-lightgrey)](https://www.rust-lang.org)
[![tests](https://img.shields.io/badge/tests-389%20passing-brightgreen)]()
[![Panaroo parity](https://img.shields.io/badge/panaroo%20parity-full-brightgreen)]()

## Features

### Core Pipeline
- **GPU-accelerated clustering** -- MMseqs2 with CUDA support, built-in CPU fallback with SIMD acceleration
- **Iterative multi-threshold collapsing** -- Progressive gene family collapsing at [0.99, 0.95, 0.9, 0.8, 0.7] thresholds (Panaroo-compatible)
- **Length-based clustering filter** -- `--len-dif-percent` rejects gene pairs with excessive length difference
- **Memory-mapped I/O** -- Zero-copy file access for large datasets
- **Concurrent graph** -- Lock-free DashMap-based graph construction with multi-centroid nodes
- **Streaming pipeline** -- Chunked bincode+zstd for datasets larger than RAM
- **Mixed input** -- Accept GFF3, FASTA (.fna/.fa), and GenBank (.gb/.gbk/.gbff) files
- **Per-genome contig-end tracking** -- HashSet-based per-genome boundary tracking replaces boolean flags

### Error Correction (7 stages, matching Panaroo's order)
- **Paralog resolution** -- Shortest-path BFS as primary method, context vectors as fallback (depth 5) + centroid length check
- **Mistranslation correction** -- 99% identity, 95% coverage edge-pair merge
- **Iterative gene family collapsing** -- 5-threshold progressive collapse [0.99, 0.95, 0.9, 0.8, 0.7] with DistanceCache reuse
- **Contig-end pruning** -- Per-genome boundary tracking removes terminal low-support nodes
- **Missing gene recovery** -- Semi-global HW alignment in flanking contig sequences + consensus removal (strict mode)
- **Re-collapse gene families** -- DistanceCache reuse from initial collapsing stage
- **Misassembly edge cleaning** -- Contig-end + disproportionate edge removal
- **Highly variable gene detection** -- Cycle-based graph analysis (Panaroo-compatible algorithm)
- **Consensus removal** -- Deletes spurious nodes where refound hits exceed original support (strict mode)

### Pre-processing & QC
- **CheckM2 integration** -- Assembly completeness and contamination scoring
- **Distance estimation** -- skani (sparse k-mer chaining, 50x faster than FastANI)
- **MDS projection** -- Pure Rust classical MDS for genome distance visualization
- **QC visualization** -- HTML report with d3.js MDS scatter + bar charts
- **Bakta re-annotation** -- Annotate raw genome assemblies before analysis
- **Prodigal gene calling** -- `--pipeline-mode prodigal` for gene prediction on raw FASTA assemblies

### Subcommands
- **Integrate mode** (`panminer integrate`) -- Add new genomes to an existing pangenome without full rebuild
- **Standalone MSA** (`panminer msa`) -- Post-run multiple sequence alignment with core/pan mode
- **Gene extraction** (`panminer extract-gene`) -- Retrieve member sequences by cluster ID
- **GFF3 output** -- Per-genome GFF3 reconstruction from corrected graph
- **Pangrowth** (`--pangrowth`) -- Exact pangenome openness estimation via Heaps' law alpha fitting
- **Abundance visualization** (`--abundance`) -- D3.js HTML report with U-shape plot, rarefaction, and partition bars

### Output Formats (15+)
- **Roary-compatible CSV** -- `gene_presence_absence.csv` (14 metadata columns)
- **Roary gene member CSV** -- `gene_presence_absence_roary.csv` (semicolon gene IDs)
- **Binary matrix** -- `gene_presence_absence.Rtab`
- **Enriched GML** -- Graph with length, seq, protein, genome_ids, member, is_paralog, is_highly_variable
- **Panaroo reference files** -- `pan_genome_reference.fa`, `gene_data.csv` (with DNA/protein + location)
- **Combined FASTA** -- `combined_DNA_CDS.fasta`, `combined_protein_CDS.fasta`
- **Core genome alignment** -- MAFFT, Clustal Omega, or PRANK
- **Alignment trimming** -- ClipKIT + BMGE filtering
- **Codon alignment** -- MACSE v2 via Java subprocess
- **Structural variant matrix** -- `struct_presence_absence.csv` and `.tsv`
- **GFF3 output** -- Per-genome GFF3 reconstruction from corrected graph
- **JSON/JSONL** -- `_pangenome.json` and `_pangenome.jsonl`
- **Parquet** -- Apache Arrow/Parquet (feature-gated)
- **HTML visualization** -- d3.js force-directed graph (feature-gated)
- **Abundance visualization** -- D3.js U-shape plot, rarefaction, partition bars (feature-gated)
- **Summary statistics** -- Core/Soft core/Shell/Cloud classification + highly variable gene count

### Downstream Analysis (`panminer analyze`)
- **Pyseer** -- Pan-GWAS via subprocess
- **Scoary2** -- Gene-trait association (Genome Biology 2024)
- **SpydrPick** -- MI-based epistasis detection (NAR 2019)
- **Panstripe** -- Phylogeny-aware gene gain/loss rates (Genome Biology 2023)
- **AMRFinderPlus** -- Curated AMR detection (NCBI)
- **Pangrowth** -- Exact pangenome openness estimation (Parmigiani, Wittler, Stoye 2024)
- **Gene neighborhood** -- Native BFS extraction from pangenome graph
- **Accumulation curves** -- Native rarefaction with Heaps' law fitting
- **GrapeTree/iTOL** -- Native profile/annotation export

### Infrastructure
- **Subprocess timeout** -- All external tools run with 1-hour timeout protection
- **Error on feature absence** -- cDBG mode returns `Error::FeatureNotEnabled` instead of silently continuing
- **cDBG pipeline mode** -- GGCAT + ggCaller for de novo gene calling
- **Prodigal mode** -- Gene prediction on raw FASTA assemblies

## Supported Platforms

**Linux** and **macOS** are the primary supported platforms. This aligns with the bioinformatics ecosystem -- MMseqs2, CD-HIT, Bakta, CheckM2, skani, and other external tools are available via conda on Linux/macOS only.

Windows users should use **WSL2** (Windows Subsystem for Linux) or run PanMiner on a Linux server/HPC cluster.

## Installation

### Prerequisites

- **Rust 1.70+** -- [Install via rustup](https://rustup.rs)
- **MMseqs2** (optional) -- For GPU-accelerated clustering (`conda install -c bioconda mmseqs2`)
- **skani** (optional) -- For fast, robust ANI distance estimation (`conda install -c bioconda skani`)
- **Bakta** (optional) -- For re-annotation of raw genome assemblies
- **Prodigal** (optional) -- For gene prediction on raw FASTA (`conda install -c bioconda prodigal`)
- **CheckM2** (optional) -- For pre-processing quality control
- **MAFFT/Clustal/PRANK** (optional) -- For multiple sequence alignment
- **ClipKIT** (optional) -- For alignment trimming
- **pyseer** (optional) -- For pan-GWAS analysis
- **Pangrowth** (optional) -- For exact pangenome openness estimation

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

The `environment.yml` installs PanMiner's external tool dependencies (MMseqs2, skani, MAFFT, etc.) but **not Rust itself** -- Rust should come from rustup to avoid conda solver conflicts.

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
conda install -c conda-forge -c bioconda mmseqs2 skani mafft clustalo prank clipkit prodigal

# Alignment filtering
pip install bmge                                  # BMGE (requires Biopython)
pip install biopython                              # Biopython dependency for BMGE

# Downstream analysis
conda install -c conda-forge -c bioconda pyseer   # Pan-GWAS
conda install -c conda-forge -c bioconda scoary2  # Gene-trait association
conda install -c conda-forge -c bioconda spydrpick # Epistasis detection
conda install -c conda-forge -c bioconda r-base=4.3 r-panstripe  # Evolutionary model
conda install -c conda-forge -c bioconda ncbi-amrfinder  # AMR detection

# Pangenome openness
conda install -c conda-forge -c bioconda pangrowth  # Exact openness estimation

# Tools with large dep trees -- install separately if needed
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

### Prodigal gene calling on raw FASTA

```bash
# Use Prodigal for gene prediction instead of Bakta
panminer *.fasta -o results --pipeline-mode prodigal
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

### Integrate new genomes into existing pangenome

```bash
panminer integrate --graph results --input new_genome.gff -o updated_results
```

### Post-run multiple sequence alignment

```bash
# Core gene alignment
panminer msa -i results --mode core --aligner mafft -o alignments/

# Full pangenome alignment
panminer msa -i results --mode pan --aligner clustal -o alignments/
```

### Extract gene sequences by cluster

```bash
# Extract DNA sequences
panminer extract-gene -i results --cluster cluster_0042 -o gene_seqs.fasta

# Extract protein sequences
panminer extract-gene -i results --cluster cluster_0042 --protein -o protein_seqs.fasta
```

### Pangenome openness estimation

```bash
panminer analyze -i results --pangrowth
```

### Abundance visualization

```bash
panminer analyze -i results --abundance
```

### Generate specific output formats

```bash
panminer *.gff -o output --formats matrix,json,graph,gff
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

### Custom iterative collapsing thresholds

```bash
# Override default 5-threshold collapsing
panminer *.gff -o output --collapse-thresholds "0.99,0.95,0.9"
```

### Length-based clustering filter

```bash
# Reject gene pairs with >5% length difference (default: 2%)
panminer *.gff -o output --len-dif-percent 0.95
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

# Pangenome openness estimation
panminer analyze -i results --pangrowth

# Abundance visualization (U-shape, rarefaction, partition bars)
panminer analyze -i results --abundance

# P/A filtering
panminer filter-pa -i results -o filtered_output
```

## Command Line Reference

### Main Pipeline (`panminer`)

| Flag | Description | Default |
|------|-------------|---------|
| `INPUT` | Input GFF3/FASTA/GenBank files (required) | -- |
| `-o, --output` | Output directory | `panminer_output` |
| `-t, --threads` | Thread count (0 = auto) | `0` |
| `--chunk-size` | Genomes per chunk for streaming | `100` |
| `--identity` | Clustering identity threshold (0.5--1.0) | `0.98` |
| `--mode` | Correction mode: strict, default, sensitive | `default` |
| `--pipeline-mode` | Pipeline mode: gff, dbg, or prodigal | `gff` |
| `--kmer-size` | K-mer size for cDBG mode (15--127) | `31` |
| `--len-dif-percent` | Length difference cutoff for clustering (0.0--1.0) | `0.98` |
| `--collapse-thresholds` | Comma-separated iterative collapsing thresholds | `0.99,0.95,0.9,0.8,0.7` |
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
| `--formats` | Output formats (comma-separated: matrix, alignment, graph, json, gff, struct, parquet, html) | `matrix,alignment,graph` |
| `--trim-alignment` | Trim alignment with ClipKIT | `false` |
| `--trim-mode` | ClipKIT trim mode | `smart-gap` |
| `--filter-alignment` | Alignment filter: none, clipkit, bmge | `none` |
| `--codons` | Generate codon alignments (MACSE) | `false` |
| `--gwas` | Run GWAS with pyseer | `false` |
| `--phenotype` | Phenotype file for GWAS | None |

### Integrate Mode (`panminer integrate`)

Add new genomes to an existing pangenome without a full rebuild.

| Flag | Description |
|------|-------------|
| `--graph <dir>` | Existing PanMiner output directory (required) |
| `--input <file>` | New genome file(s) to add (required) |
| `-o, --output <dir>` | Output directory for updated pangenome |
| `--identity` | Clustering identity threshold for matching |
| `-t, --threads` | Thread count (0 = auto) |
| `--verbose` | Enable debug logging |

```bash
panminer integrate --graph results --input new_genome.gff -o updated_results
panminer integrate --graph results --input *.gff -o updated_results --identity 0.95
```

### Gene Extraction (`panminer extract-gene`)

Retrieve member sequences from a pangenome by cluster ID.

| Flag | Description |
|------|-------------|
| `-i, --input <dir>` | PanMiner output directory (required) |
| `--cluster <id>` | Cluster ID to extract (required) |
| `-o, --output <file>` | Output FASTA file |
| `--protein` | Extract protein sequences instead of DNA |

```bash
panminer extract-gene -i results --cluster cluster_0042 -o gene_seqs.fasta
panminer extract-gene -i results --cluster cluster_0042 --protein -o protein_seqs.fasta
```

### Standalone MSA (`panminer msa`)

Post-run multiple sequence alignment with core or pan mode.

| Flag | Description |
|------|-------------|
| `-i, --input <dir>` | PanMiner output directory (required) |
| `-o, --output <dir>` | Output directory for alignments |
| `--mode` | Alignment mode: core or pan (default: core) |
| `--aligner` | Aligner: mafft, clustal, or prank (default: mafft) |
| `-t, --threads` | Thread count (0 = auto) |
| `--verbose` | Enable debug logging |

```bash
panminer msa -i results --mode core --aligner mafft -o alignments/
panminer msa -i results --mode pan --aligner clustal -o alignments/
```

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
| `--pangrowth` | Run Pangrowth openness estimation |
| `--abundance` | Generate abundance visualization |

### P/A Filtering (`panminer filter-pa`)

| Flag | Description |
|------|-------------|
| `-i, --input <dir>` | PanMiner output directory (required) |
| `-o, --output <dir>` | Filtered output directory |
| `--min-genomes <n>` | Minimum number of genomes a gene must be in |
| `--max-genomes <n>` | Maximum number of genomes a gene must be in |
| `--remove-fragments` | Remove fragmented genes |
| `--remove-pseudogenes` | Remove pseudogenes |

## Algorithm Details

### Correction Pipeline (7 stages, matching Panaroo's order)

The correction pipeline follows Panaroo's exact algorithmic order to ensure compatibility:

| Stage | Algorithm | Description |
|-------|-----------|-------------|
| 4.1 | Paralog resolution | Shortest-path BFS as primary method, context vector similarity (depth 5) as fallback + centroid length check |
| 4.2 | Mistranslation correction | Merge adjacent clusters where coverage >= 95% and identity >= 99% |
| 4.3 | Iterative gene family collapsing | Progressive collapse at [0.99, 0.95, 0.9, 0.8, 0.7] thresholds; DistanceCache reused across thresholds |
| 4.4 | Contig-end pruning | Per-genome HashSet boundary tracking; removes terminal low-support nodes |
| 4.5 | Missing gene recovery | Semi-global HW alignment in 5000bp flanking regions; consensus removal in strict mode deletes spurious nodes where refound hits exceed original support |
| 4.6 | Re-collapse gene families | DistanceCache reuse from stage 4.3; cleans up graph after gene recovery |
| 4.7 | Misassembly edge cleaning | Two-criteria removal: contig-end edges + disproportionate edge support |

### Correction Modes

| Parameter | Strict | Default | Sensitive |
|-----------|--------|---------|-----------|
| Contamination threshold | max(2, ceil(0.05N)) | 2 | 1 |
| Consensus removal | Yes | No | No |
| Edge support | max(2, ceil(0.01N)) | max(2, ceil(0.01N)) | 0 (disabled) |
| Trailing recursion | Unlimited | Unlimited | Disabled |
| Use case | Phylogenetic studies | General analysis | Preserve rare genes/plasmids |

N = number of genomes. Strict mode aggressively removes contamination and spurious refound genes, making it best for phylogenetic studies. Sensitive mode preserves rare elements like plasmids and low-abundance genes.

### Clustering

**MMseqs2 GPU:**
- Uses `easy-cluster` with `--cov-mode 1 -c {len_dif_percent}` for length-based filtering
- Auto-detects GPU via `mmseqs version` output (searches for CUDA/GPU strings)
- Passes `--gpu 1` when GPU is available

**CPU fallback:**
- Greedy incremental clustering with SIMD-accelerated sequence comparison (AVX2/NEON/scalar with runtime detection)
- Length difference filtering: rejects pairs where length ratio falls below `--len-dif-percent` (default: 0.98)

Both clusterers support `--len-dif-percent` to prevent merging genes with excessive length differences, reducing false merges of truncated or fragmented sequences.

### Pangenome Openness

PanMiner uses **Pangrowth** (Parmigiani, Wittler, Stoye 2024) for exact pangenome growth and core curves:

- Computes exact pangenome size and core genome size for all genome subsets
- Heaps' law fitting: pangenome size = k x n^alpha
- alpha > 0 indicates an **open** pangenome (new genes continue to appear)
- alpha <= 0 indicates a **closed** pangenome (finite gene pool)
- Does not require a phylogenetic tree (unlike Panstripe)

## Correction Modes

| Mode | Contamination Threshold | Use Case |
|------|------------------------|----------|
| `strict` | 5% of genome count | Phylogenetic studies -- removes more, enables consensus removal |
| `default` | 2 genomes | Most use cases -- balanced |
| `sensitive` | 1 genome | Preserve rare genes/plasmids -- minimal correction |

## Distance Estimation

PanMiner uses **skani** for distance estimation -- the fastest and most robust ANI tool available.

1. **skani** -- 50x faster than FastANI, robust to incomplete genomes (MAGs)

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
| `*.gff` | Per-genome GFF3 reconstruction from corrected graph | No |

### QC Output Files (when enabled)

| File | Description |
|------|-------------|
| `qc_stats.csv` | Per-genome QC metrics (completeness, contamination, ANI, etc.) |
| `qc_summary.txt` | Human-readable QC summary with pass/fail status |
| `qc_viz.html` | Interactive HTML report with MDS scatter + bar charts |

### Visualization & Analysis Files

| File | Description |
|------|-------------|
| `_pangenome.json` | Full pangenome summary (JSON) |
| `_pangenome.jsonl` | Streaming JSON output (one cluster per line) |
| `_pangenome.parquet` | Parquet format (requires `--features parquet`) |
| `_pangenome.html` | Interactive HTML force-directed graph (requires `--features viz`) |
| `abundance_viz.html` | Abundance visualization with U-shape plot, rarefaction, partition bars (requires `--features viz`) |

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
Phase 0:    QC (optional) -- CheckM2 completeness/contamination + skani ANI distance
Phase 0.5:  Re-annotation (optional) -- Bakta or Prodigal gene calling
Phase 1:    Parse -- Memory-mapped GFF3/FASTA/GenBank parsing with Rayon parallelism
Phase 2:    Cluster -- MMseqs2-GPU or CPU with length-based filtering (--len-dif-percent)
Phase 3:    Build graph -- DashMap concurrent graph (multi-centroid nodes + per-genome contig-end tracking)
Phase 4:    Correct -- 7-stage correction pipeline:
             4.1 Paralog resolution (shortest-path BFS + context vectors)
             4.2 Mistranslation correction (99% identity, 95% coverage)
             4.3 Iterative gene family collapsing [0.99, 0.95, 0.9, 0.8, 0.7]
             4.4 Contig-end pruning (per-genome tracking)
             4.5 Missing gene recovery + consensus removal (strict mode)
             4.6 Re-collapse gene families (DistanceCache reuse)
             4.7 Misassembly edge cleaning
Phase 4.8:  Highly variable gene detection (cycle-based)
Phase 5:    Matrix -- Build presence/absence matrix (BitPackedMatrix, 8x memory reduction)
Phase 6:    Output -- Generate 15+ output formats
Phase 7:    Downstream (optional) -- GWAS, evolution, AMR, visualization
```

## Architecture

```
src/
├── lib.rs              # Library entry point, re-exports
├── main.rs             # CLI with clap (main + integrate + extract-gene + msa subcommands)
├── config.rs           # Configuration structs (CorrectionMode, OutputFormat, PipelineMode, FilterMethod)
├── error.rs            # Error types (thiserror-based)
├── pipeline.rs         # Main pipeline orchestration (PanminerPipeline + cDBG mode)
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
│   ├── ggcaller.rs     # ggCaller gene calling subprocess
│   ├── orphos.rs       # Prodigal gene calling (feature-gated)
│   ├── extract_gene.rs # Gene extraction subcommand
│   └── integrate.rs    # Integrate mode
│
├── clustering/         # Compute -- Clustering
│   ├── traits.rs       # Clusterer trait
│   ├── mmseqs.rs       # MMseqs2-GPU integration (with --len-dif-percent)
│   ├── cpu.rs          # CPU fallback with SIMD + length filtering
│   └── alignment_*.rs  # MSA tools (MAFFT, Clustal, PRANK)
│
├── graph/              # Data Structures
│   ├── types.rs        # Gene, GeneCluster, Node (multi-centroid), Edge, PangenomeGraph
│   ├── concurrent.rs   # DashMap-based ConcurrentGraph
│   ├── matrix.rs       # BitPackedMatrix (8x memory reduction)
│   ├── builder.rs      # Graph construction with gene_members + gene_lookup
│   ├── structural_variants.rs  # Structural variant detection
│   ├── highly_variable.rs     # Highly variable gene detection
│   └── merge.rs        # Pangenome merging
│
├── correction/         # Error Correction (7 stages)
│   ├── contamination.rs    # Low-support node removal
│   ├── contig_end.rs       # Contig-end pruning (per-genome HashSet tracking)
│   ├── fragment.rs        # Fragment merging + iterative collapsing with DistanceCache
│   ├── missing.rs          # Missing gene recovery (semi-global HW) + consensus removal
│   ├── misassembly.rs     # Misassembly edge cleaning
│   ├── paralog.rs          # Paralog resolution (shortest-path BFS + context vectors)
│   └── simd.rs             # SIMD sequence comparison (AVX2/NEON) + Levenshtein
│
├── downstream/         # Downstream Analysis
│   ├── gwas/            # PyseerRunner, Scoary2Runner, SpydrPickRunner
│   ├── evolution/       # PanstripeRunner, PangrowthRunner
│   ├── resistome/       # AmrFinderRunner
│   └── exploration/     # AccumulationCurve, GeneNeighborhood, GrapeTree, iTOL
│
├── gwas/               # Top-level GWAS module
│   ├── traits.rs       # GWASRunner trait
│   └── pyseer.rs        # PyseerRunner
│
└── output/              # Output Generation (17+ files)
    ├── mod.rs           # OutputWriter + OutputPaths
    ├── matrix.rs        # Roary CSV (14-col), Rtab, gene member CSV
    ├── alignment.rs     # Core/accessory alignment via MSA
    ├── graph.rs         # GML (enriched attributes)
    ├── json.rs          # gene_data.csv, pan_genome_reference.fa, JSON
    ├── gff.rs           # Per-genome GFF3 reconstruction
    ├── struct_csv.rs    # Structural variant CSV
    ├── sv_matrix.rs     # Structural variant TSV
    ├── summary.rs       # Core/Soft core/Shell/Cloud + highly variable
    ├── parquet.rs       # Apache Parquet (feature-gated)
    ├── html_viz.rs      # d3.js HTML visualization (feature-gated)
    ├── abundance_viz.rs # Abundance visualization (feature-gated)
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
| Iterative collapsing | 2-pass (0.98, 0.7) | 5-threshold [0.99--0.7] |
| Length filtering | No | Yes (--len-dif-percent) |
| Graph | NetworkX | DashMap concurrent + petgraph |
| Distance | Mash | skani (sparse k-mer chaining) |
| I/O | Standard reads | Memory-mapped (mmap2) |
| Memory | High (NetworkX) | Low (BitPackedMatrix 8x reduction) |
| Streaming | No | Yes (chunked bincode+zstd) |
| GPU | None | MMseqs2-GPU subprocess |
| Parallelism | multiprocessing | Rayon work-stealing |
| Pipeline modes | GFF3 only | GFF3 + cDBG + Prodigal |
| Error correction | 6 modules | 7 stages + highly variable detection |
| Consensus removal | Yes (strict) | Yes (strict) |
| Shortest-path paralog | Yes (nx.shortest_path) | Yes (petgraph BFS) |
| Multi-centroid nodes | Yes (centroid: list) | Yes (Vec<Sequence>) |
| Per-genome contig-end | Yes (hasEnd per gene) | Yes (HashSet<GenomeId>) |
| Downstream tools | Scoary, SpydrPick | Scoary2, SpydrPick, Panstripe, AMRFinderPlus |
| Integrate mode | No | Yes (panminer integrate) |
| Standalone MSA | Yes (panaroo-msa) | Yes (panminer msa) |
| Gene extraction | Yes (extract-gene) | Yes (panminer extract-gene) |
| GFF output | No | Yes (--formats gff) |
| Pangenome openness | No | Yes (Pangrowth) |
| Abundance viz | No | Yes (--abundance) |
| Prodigal mode | No | Yes (--pipeline-mode prodigal) |

## Build Features

| Feature | Description |
|---------|-------------|
| `cpu` (default) | CPU processing support |
| `mmseqs` | MMseqs2 integration for clustering |
| `parquet` | Parquet output format support |
| `python` | PyO3 Python bindings |
| `viz` | HTML visualization output |
| `dbg` | GGCAT cDBG pipeline mode |
| `prodigal` | Prodigal gene calling support (subprocess) |
| `full` | All features enabled (cpu, mmseqs, parquet, viz, dbg, prodigal) |

```bash
# Build with specific features
cargo build --features mmseqs,parquet

# Build with all features
cargo build --features full

# Build with Prodigal support
cargo build --features prodigal
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

MIT License -- see [LICENSE](LICENSE) for details.

## Citation

If you use PanMiner in your research, please cite:

- PanMiner: [DOI/Preprint link TBD]
- Panaroo: Tonkin-Hill et al. (2020) "Rethinking the pangenome: high-quality bacterial genomes made easy using Panaroo" *Genome Biology* 21:258
- MMseqs2: Steinegger & Soding (2017) "MMseqs2 enables sensitive protein sequence searching for the massive scale of sequence data" *Nature Biotechnology* 35:1026
- skani: Shaw & Yu (2023) "skani: fast, lightweight nucleotide comparison with accurate pairwise ANI estimation" *Nature Methods* 20:1160
- Pangrowth: Parmigiani, Wittler, Stoye (2024) "Pangrowth: The dynamics of pangenomes" *PCI Computational Biology*

## Acknowledgments

- PanMiner is inspired by [Panaroo](https://github.com/gtonkinhill/panaroo) (Tonkin-Hill et al., Genome Biology 2020)
- Built with Rust, Rayon, DashMap, memmap2, and petgraph
- Uses MMseqs2 for GPU-accelerated clustering (Steinegger & Soding, Nature Biotechnology 2017)
- Uses skani for fast ANI estimation (Shaw & Yu, Nature Methods 2023)
- Uses Bakta for genome re-annotation (Schwengers et al., Microbial Genomics 2021)
- Uses Pangrowth for exact pangenome openness estimation (Parmigiani, Wittler, Stoye, PCI Computational Biology 2024)