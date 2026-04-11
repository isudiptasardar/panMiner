# Downstream Analysis Module Design

**Date**: 2026-04-12
**Status**: Draft
**Scope**: Trait-based downstream analysis module with subprocess wrappers for GWAS, evolutionary modeling, AMR detection, and pangenome exploration

---

## 1. Overview

PanMiner currently has strong pangenome construction and error correction but lacks downstream analysis capabilities that Panaroo and modern pangenome tools provide. This design adds a composable downstream analysis module following PanMiner's existing trait-based architecture pattern.

### Goals

- Match and exceed Panaroo's downstream analysis subcommands
- Integrate the best modern tools (Panstripe over IMG/FMG, Scoary2 over Scoary, AMRFinderPlus)
- Maintain PanMiner's trait-based architecture for extensibility
- Use subprocess wrappers (matching existing pattern for pyseer, MAFFT, ClipKIT)
- Add native Rust implementations where no external tool is needed (neighborhood extraction, accumulation curves)

### Non-Goals

- Reimplementing algorithms natively in Rust (use subprocess wrappers instead)
- Building a workflow engine (users compose analyses via CLI flags)
- Supporting eukaryotic pangenomes (focus on prokaryotic use cases)

---

## 2. Architecture

### 2.1 Module Structure

```
src/downstream/
├── mod.rs              # Re-exports, DownstreamRunner trait
├── traits.rs           # DownstreamRunner trait + DownstreamResult
├── gwas/
│   ├── mod.rs          # Re-exports
│   ├── pyseer.rs       # PyseerRunner (moved from src/gwas/)
│   ├── scoary.rs       # Scoary2Runner — gene-trait association
│   └── spydrpick.rs    # SpydrPickRunner — co-selection/epistasis
├── evolution/
│   ├── mod.rs          # Re-exports
│   └── panstripe.rs    # PanstripeRunner — GLM gene gain/loss rates
├── resistome/
│   ├── mod.rs          # Re-exports
│   └── amrfinder.rs    # AmrFinderRunner — AMR gene detection
└── exploration/
    ├── mod.rs          # Re-exports
    ├── neighborhood.rs # GeneNeighborhoodExtractor — graph neighborhood (native)
    ├── accumulation.rs # AccumulationCurveRunner — gene accumulation + Heaps' law (native)
    └── grapetree.rs    # GrapeTreeExportRunner — export for GrapeTree/iTOL
```

### 2.2 Core Trait

```rust
/// Trait for all downstream analysis runners.
pub trait DownstreamRunner: Send + Sync {
    type Output: DownstreamResult;

    /// Run the analysis given an output directory containing PanMiner outputs.
    /// Runners read what they need from the output directory (graph GML, P/A CSV,
    /// protein FASTA, etc.) rather than requiring in-memory graph+matrix.
    fn run(&self, output_dir: &Path) -> Result<Self::Output>;

    /// Name of the tool/analysis.
    fn name(&self) -> &str;

    /// Check if the external tool is installed and available.
    fn is_available(&self) -> bool;

    /// Declare required input files that must exist in the output directory.
    fn required_inputs(&self) -> Vec<DownstreamInput>;
}

/// Input types that downstream analyses may require from the output directory.
pub enum DownstreamInput {
    FinalGraph,           // final_graph.gml
    PresenceAbsenceCsv,   // gene_presence_absence.csv
    ProteinFasta,         // combined_protein_CDS.fasta
    DnaFasta,             // combined_DNA_CDS.fasta
    GeneDataCsv,          // gene_data.csv
    PhenotypesFile,       // user-supplied phenotypes
    PhylogeneticTree,     // user-supplied Newick tree
    AmrDatabase,          // AMRFinderPlus database path
}

/// Result trait for downstream analysis outputs.
pub trait DownstreamResult: Send + Sync {
    /// Write all output files to the specified directory.
    fn write_to(&self, dir: &Path) -> Result<()>;

    /// Summary statistics as a string.
    fn summary(&self) -> String;
}
```

### 2.3 CLI Integration

Add a new `panminer analyze` subcommand to `main.rs`:

```rust
// In main.rs, add to the clap App:
.subcommand(
    App::new("analyze")
        .about("Run downstream analyses on pangenome output")
        .arg(Arg::new("input").short('i').long("input").required(true))
        .arg(Arg::new("gwas").long("gwas").help("Run GWAS analysis"))
        .arg(Arg::new("gwas-tool").long("gwas-tool").default_value("pyseer")
             .possible_values(["pyseer", "scoary2", "spydrpick"]))
        .arg(Arg::new("phenotypes").long("phenotypes").help("Path to phenotypes file"))
        .arg(Arg::new("panstripe").long("panstripe").help("Gene gain/loss rate estimation"))
        .arg(Arg::new("tree").long("tree").help("Path to phylogenetic tree (Newick)"))
        .arg(Arg::new("amr").long("amr").help("AMR gene detection"))
        .arg(Arg::new("amr-database").long("amr-database").help("AMRFinderPlus database path"))
        .arg(Arg::new("neighborhood").long("neighborhood").help("Gene neighborhood extraction"))
        .arg(Arg::new("seed-gene").long("seed-gene").help("Seed gene/cluster ID"))
        .arg(Arg::new("neighborhood-depth").long("neighborhood-depth").default_value("5"))
        .arg(Arg::new("accumulation").long("accumulation").help("Gene accumulation curves"))
        .arg(Arg::new("num-samples").long("num-samples").default_value("100"))
        .arg(Arg::new("export-grapetree").long("export-grapetree").help("Export for GrapeTree"))
        .arg(Arg::new("export-itol").long("export-itol").help("Export for iTOL"))
)
```

---

## 3. Runner Specifications

### 3.1 GWAS Module

#### 3.1.1 PyseerRunner (relocated from `src/gwas/`)

**Status**: Already implemented. Needs relocation and enhancement.

**Changes from current implementation**:
- Move from `src/gwas/pyseer.rs` to `src/downstream/gwas/pyseer.rs`
- Implement `DownstreamRunner` trait (the existing `GWASRunner` trait in `src/gwas/traits.rs` will be kept for backward compatibility; `DownstreamRunner` will delegate to `GWASRunner` internally)
- Add `--unitigs` input mode support
- Improve phenotype file generation (currently hardcoded to `total_genes`)
- Add `--lmm` (linear mixed model) flag passthrough
- Add `--burden` testing mode passthrough

**External dependency**: `pyseer` (Python, install via pip)

#### 3.1.2 Scoary2Runner

**Status**: New implementation.

**Purpose**: Gene-trait association testing. Scoary2 extends Scoary with support for continuous phenotypes, multi-omics data, and an interactive data exploration app.

**Implementation**:
```rust
pub struct Scoary2Runner {
    scoary2_path: PathBuf,
    phenotypes_file: Option<PathBuf>,
    output_dir: Option<PathBuf>,
    threads: usize,
}
```

**Algorithm**:
1. Write P/A matrix to temporary CSV in Scoary2 input format
2. Write phenotypes file (binary or continuous)
3. Run `scoary2 -t <phenotypes> -p <pa_matrix> -o <output_dir>`
4. Parse output CSV: gene IDs, trait associations, p-values, FDR-corrected p-values
5. Write results to `scoary_results.csv`

**External dependency**: `scoary2` (Python, install via pip)

**Key reference**: Roder et al. (2024) *Genome Biology* 25:93

#### 3.1.3 SpydrPickRunner

**Status**: New implementation.

**Purpose**: Identify gene presence/absence patterns that are highly correlated or anti-correlated while accounting for population structure.

**Implementation**:
```rust
pub struct SpydrPickRunner {
    spydrpick_path: PathBuf,
    output_dir: Option<PathBuf>,
}
```

**Algorithm**:
1. Write P/A matrix to temporary TSV
2. Run `spydrpick -i <pa_matrix> -o <output_dir>`
3. Parse output: gene pair correlations, MI scores, p-values
4. Write results to `spydrpick_correlations.csv`

**External dependency**: `spydrpick` (C++, install via conda or compile)

**Key reference**: Pensar et al. (2019) *Nucleic Acids Research*

---

### 3.2 Evolution Module

#### 3.2.1 PanstripeRunner

**Status**: New implementation.

**Purpose**: Estimate gene gain and loss rates using phylogenetically informed generalized linear models. **This replaces Panaroo's IMG/FMG models** — Panstripe (2023) outperforms both in simulations, especially under annotation errors and sampling bias.

**Why Panstripe over IMG/FMG**:
- Robust to annotation/clustering errors (critical for real data)
- Controls for population structure and sampling bias
- Includes a "tip" covariate to absorb spurious gene calls at terminal branches
- Can compare rates between pangenomes and test associations with covariates
- Published benchmarking shows superior performance

**Implementation**:
```rust
pub struct PanstripeRunner {
    rscript_path: PathBuf,
    panstripe_lib_path: Option<PathBuf>,
    tree_file: Option<PathBuf>,  // Newick phylogenetic tree
    output_dir: Option<PathBuf>,
}
```

**Algorithm**:
1. Detect R and panstripe R package
2. Write P/A matrix and phylogenetic tree to temporary files
3. Run R script: `panstripe::run_panstripe(pa_matrix, tree, ...)`
4. Parse results: gain rate, loss rate, pangenome openness (alpha), diagnostic plots
5. Write `panstripe_rates.csv` and `panstripe_diagnostics.pdf`

**External dependency**: R + `panstripe` package (install via `install.packages("panstripe")`)

**Key reference**: Tonkin-Hill et al. (2023) *Genome Research* 33(1):129–140

---

### 3.3 Resistome Module

#### 3.3.1 AmrFinderRunner

**Status**: New implementation.

**Purpose**: Detect antimicrobial resistance genes, stress response genes, and virulence factors in genome assemblies.

**Why AMRFinderPlus over CARD/RGI**:
- NCBI-curated with continuously updated database
- Public domain (no commercial license needed)
- Hierarchical detection: EXACT > ALLELE > BLAST > HMM > PARTIAL > POINT > INTERNAL_STOP
- Taxon-specific point mutations
- Evidence tracking (method, identity, coverage, contig position)
- Integrated into NCBI Pathogen Detection pipeline (validated on 800K+ isolates)

**Implementation**:
```rust
pub struct AmrFinderRunner {
    amrfinder_path: PathBuf,
    database_path: Option<PathBuf>,
    organism: Option<String>,    // Taxon-specific analysis
    threads: usize,
    output_dir: Option<PathBuf>,
}
```

**Algorithm**:
1. Detect `amrfinder` binary via `which`
2. Write protein FASTA (from `combined_protein_CDS.fasta`) and GFF annotations to temp files
3. Run `amrfinder_plus -i <input_dir> -o <output> --plus` (Plus mode includes stress/virulence)
4. Optionally add `--organism <species>` for taxon-specific point mutations
5. Parse TSV output: gene name, drug class, method, identity, coverage, etc.
6. Write `amr_results.tsv` and `amr_summary.txt`

**External dependency**: `amrfinder` (NCBI AMRFinderPlus, install via conda or NCBI FTP)

**Key reference**: Feldgarden et al. (2021) *Scientific Reports* 11:12720

---

### 3.4 Exploration Module

#### 3.4.1 GeneNeighborhoodExtractor (Native Rust)

**Status**: New implementation. **No external dependency** — operates directly on PangenomeGraph.

**Purpose**: Extract genes within N hops of a seed gene in the pangenome graph. Mirrors `panaroo-gene-neighbourhood`.

**Implementation**:
```rust
pub struct GeneNeighborhoodExtractor {
    seed_gene: ClusterId,
    max_depth: usize,  // default: 5
    output_dir: Option<PathBuf>,
}
```

**Algorithm**:
1. Find the seed gene cluster in PangenomeGraph
2. BFS from seed node, collecting all nodes within `max_depth` hops
3. Record: cluster IDs, support, annotations, genome membership, hop distance
4. Extract subgraph edges
5. Write `neighborhood_genes.csv` (gene info + hop distance)
6. Write `neighborhood_subgraph.gml` (Cytoscape-compatible)

#### 3.4.2 AccumulationCurveRunner (Native Rust)

**Status**: New implementation. **No external dependency** — operates directly on BitPackedMatrix.

**Purpose**: Generate gene accumulation curves and fit Heaps' law to classify pangenomes as open or closed. Mirrors `panaroo-plot-abundance`.

**Implementation**:
```rust
pub struct AccumulationCurveRunner {
    num_samples: usize,       // default: 100
    rarefaction_points: usize, // default: 20
    output_dir: Option<PathBuf>,
}
```

**Algorithm**:
1. Random subsampling of genomes (without replacement) at increasing genome counts
2. For each sample size k (1, 2, ..., n genomes), count total, core, and accessory genes
3. Average across `num_samples` iterations per k
4. Fit Heaps' law: n(k) = A * k^alpha using least squares on log-transformed data
5. Classify: alpha < 1.0 → closed pangenome, alpha >= 1.0 → open pangenome
6. Write `accumulation_curve.csv` (genome_count, total_genes, core_genes, accessory_genes, mean, stderr)
7. Write `heaps_law_fit.csv` (alpha, A, r_squared, openness_classification)

#### 3.4.3 GrapeTreeExportRunner

**Status**: New implementation. Subprocess wrapper for GrapeTree + native iTOL export.

**Purpose**: Export pangenome data for visualization in GrapeTree (minimum spanning trees) and iTOL (phylogenetic trees with annotation tracks).

**Implementation**:
```rust
pub struct GrapeTreeExportRunner {
    grapetree_path: Option<PathBuf>,  // optional, for MSTree computation
    output_dir: Option<PathBuf>,
}
```

**Algorithm**:
1. Convert P/A matrix to allelic profile format (GrapeTree input)
2. If GrapeTree installed: run `grapetree -i <profiles> -o <output>`
3. Generate iTOL annotation file from genome metadata and cluster data
4. Write `grapetree_profiles.tsv` (allelic profiles)
5. Write `itol_annotations.txt` (iTOL dataset file with colorstrip + text labels)

**External dependency**: `grapetree` (Python, optional — profiles exported regardless)

---

## 4. Output File Specification

All downstream analysis outputs are written to `<output_dir>/downstream/`.

| File | Source | Format | Description |
|------|--------|--------|-------------|
| `pyseer_associations.csv` | PyseerRunner | CSV | Gene-phenotype associations with p-values, FDR, effect sizes |
| `scoary_results.csv` | Scoary2Runner | CSV | Gene-trait associations with Fisher's exact test p-values |
| `spydrpick_correlations.csv` | SpydrPickRunner | CSV | Gene co-selection pairs with MI scores |
| `panstripe_rates.csv` | PanstripeRunner | CSV | Gene gain rate, loss rate, pangenome openness statistics |
| `panstripe_diagnostics.pdf` | PanstripeRunner | PDF | Diagnostic plots (accumulation curve, gain/loss rate estimates) |
| `amr_results.tsv` | AmrFinderRunner | TSV | AMR genes, drug classes, evidence types, identity, coverage |
| `amr_summary.txt` | AmrFinderRunner | TXT | Human-readable AMR summary by drug class |
| `neighborhood_genes.csv` | NeighborhoodExtractor | CSV | Genes in N-hop neighborhood with hop distance |
| `neighborhood_subgraph.gml` | NeighborhoodExtractor | GML | Subgraph for Cytoscape |
| `accumulation_curve.csv` | AccumulationCurve | CSV | Rarefaction curve data (genome count vs gene count) |
| `heaps_law_fit.csv` | AccumulationCurve | CSV | Alpha parameter, A coefficient, r-squared, openness |
| `grapetree_profiles.tsv` | GrapeTreeExport | TSV | Allelic profiles for GrapeTree |
| `itol_annotations.txt` | iTOLExport | TXT | iTOL dataset file with colorstrip + text labels |

---

## 5. Dependency Management

Each runner follows PanMiner's existing pattern:

```rust
impl DownstreamRunner for Scoary2Runner {
    fn is_available(&self) -> bool {
        Self::detect().is_some()
    }
}

impl Scoary2Runner {
    pub fn detect() -> Option<Self> {
        which::which("scoary2").ok().map(|path| Self::new(path))
    }
}
```

**External tool detection summary**:

| Tool | Detection | Install Method |
|------|-----------|----------------|
| `pyseer` | `which::which("pyseer")` | `pip install pyseer` |
| `scoary2` | `which::which("scoary2")` | `pip install scoary2` |
| `spydrpick` | `which::which("spydrpick")` | conda or compile from source |
| R + `panstripe` | `which::which("Rscript")` + R package check | `install.packages("panstripe")` |
| `amrfinder` | `which::which("amrfinder")` | conda or NCBI FTP |
| `grapetree` | `which::which("grapetree")` | `pip install grapetree` |

---

## 6. Implementation Priority

### Phase 1 (P1): Core GWAS + Exploration
1. Relocate PyseerRunner to downstream/gwas/ and implement DownstreamRunner trait
2. Scoary2Runner
3. GeneNeighborhoodExtractor (native)
4. AccumulationCurveRunner (native)
5. `panminer analyze` CLI subcommand

### Phase 2 (P2): Evolution + Resistome
6. PanstripeRunner
7. AmrFinderRunner
8. SpydrPickRunner
9. GrapeTreeExportRunner

### Phase 3 (P3): Enhancements
10. Pyseer unitig mode support
11. iTOL phylogenetic tree annotation
12. Accumulation curve visualization (native SVG)
13. AMR summary HTML report

---

## 7. Testing Strategy

Each runner requires two test levels:

**Unit tests** (no external tool required):
- Input file generation (P/A matrix, phenotypes, profiles)
- Output parsing logic
- Accumulation curve math (Heaps' law fitting)
- BFS neighborhood extraction
- GrapeTree/iTOL format generation

**Integration tests** (require external tool):
- Skip if tool not installed (`if !Scoary2Runner.detect() { return; }`)
- Test full pipeline: generate input → run tool → parse output → verify

---

## 8. Key Differentiators vs Panaroo

| Feature | Panaroo | PanMiner (after this design) |
|---------|---------|-------------------------------|
| GWAS tools | pyseer + Scoary | pyseer + **Scoary2** (continuous phenotypes) + **SpydrPick** |
| Evolutionary models | IMG + FMG | **Panstripe** (GLM, superior error robustness) |
| AMR detection | None | **AMRFinderPlus** (NCBI-curated, evidence-tracked) |
| Gene neighborhood | Built-in | Built-in (BFS on graph) |
| Accumulation curves | Built-in | Built-in (Heaps' law fitting) |
| Codon alignment | None | **MACSE v2** (frameshift-aware) |
| Alignment trimming | BMGE | **ClipKIT** (better benchmarking) |
| Visualization | Cytoscape GML | **D3.js HTML + GrapeTree + iTOL export** |
| Progressive updates | `panaroo-integrate` | **Graph merging** (already implemented) |