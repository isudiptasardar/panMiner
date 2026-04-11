# Downstream Analysis Implementation Plan

**Date**: 2026-04-12
**Spec**: `docs/superpowers/specs/2026-04-12-downstream-analysis-design.md`
**Phases**: 3 (independent tasks run in parallel within each phase)

---

## Phase 1: Core Infrastructure + Native Exploration Tools

### Task 1.1: Create `src/downstream/` module skeleton

**What**: Create the module directory structure, `mod.rs`, `traits.rs`, and four sub-module `mod.rs` files.

**Files to create**:
- `src/downstream/mod.rs` — Re-exports all downstream runners, declares `pub mod gwas`, `pub mod evolution`, `pub mod resistome`, `pub mod exploration`
- `src/downstream/traits.rs` — `DownstreamRunner` trait, `DownstreamResult` trait, `DownstreamInput` enum
- `src/downstream/gwas/mod.rs` — Empty re-export module (placeholder until Task 2.1)
- `src/downstream/evolution/mod.rs` — Empty re-export module (placeholder until Task 2.6)
- `src/downstream/resistome/mod.rs` — Empty re-export module (placeholder until Task 2.7)
- `src/downstream/exploration/mod.rs` — Re-exports `GeneNeighborhoodExtractor`, `AccumulationCurveRunner`, `GrapeTreeExportRunner`

**Key types** (from traits.rs):
```rust
// traits.rs
pub enum DownstreamInput {
    FinalGraph, PresenceAbsenceCsv, ProteinFasta, DnaFasta,
    GeneDataCsv, PhenotypesFile, PhylogeneticTree, AmrDatabase,
}

pub trait DownstreamRunner: Send + Sync {
    type Output: DownstreamResult;
    fn run(&self, output_dir: &Path) -> Result<Self::Output>;
    fn name(&self) -> &str;
    fn is_available(&self) -> bool;
    fn required_inputs(&self) -> Vec<DownstreamInput>;
}

pub trait DownstreamResult: Send + Sync {
    fn write_to(&self, dir: &Path) -> Result<()>;
    fn summary(&self) -> String;
}
```

**Verification**: `cargo check` compiles the new module without errors.

---

### Task 1.2: Implement `GeneNeighborhoodExtractor` (native)

**What**: BFS-based gene neighborhood extraction from `PangenomeGraph`. No external dependency.

**File**: `src/downstream/exploration/neighborhood.rs` (new file)

**Algorithm**:
1. Parse GML from `output_dir / "final_graph.gml"` → `PangenomeGraph`
2. Find the seed cluster node by `cluster_id`
3. BFS from seed node to `max_depth` hops, recording hop distance per reachable node
4. Collect subgraph edges among visited nodes
5. Write `neighborhood_genes.csv` (cluster_id, support, annotation, hop_distance, num_genomes)
6. Write `neighborhood_subgraph.gml` (Cytoscape-compatible, with hop_distance as node attribute)

**Implementation notes**:
- Use `petgraph::Graph` + BFS via `VecDeque` for neighbor traversal
- GML parsing: follow the simple line-by-line pattern from `src/graph/merge.rs::load_gml_graph()`
- Hop distance recorded as `distance` node attribute in GML output
- `NeighborhoodResult` struct implementing `DownstreamResult`

**Verification**: Unit test with synthetic graph — BFS depth 3 finds expected neighbors.

---

### Task 1.3: Implement `AccumulationCurveRunner` (native)

**What**: Rarefaction-based gene accumulation curve with Heaps' law fitting. Operates on `BitPackedMatrix`. No external dependency.

**File**: `src/downstream/exploration/accumulation.rs` (new file)

**Algorithm**:
1. Parse P/A matrix from `output_dir / "gene_presence_absence.csv"` or read from `BitPackedMatrix` if passed in-memory
2. For each k in `rare_fraction_points` evenly spaced genome counts (1 to n):
   - Randomly sample k genomes (without replacement), repeat `num_samples` times
   - Count total genes, core genes (present in all k), accessory genes
   - Compute mean and standard error across samples
3. Fit Heaps' law: `n(k) = A * k^alpha` via linear regression on log-transformed data
4. Classify: alpha >= 1.0 → "open pangenome", alpha < 1.0 → "closed pangenome"
5. Write `accumulation_curve.csv` (k, mean_total, mean_core, mean_accessory, stderr)
6. Write `heaps_law_fit.csv` (alpha, A, r_squared, classification)

**Implementation notes**:
- Use `BitPackedMatrix` methods: `count_present()`, `num_genomes()`, `num_clusters()`, `is_core()`
- Random sampling with a seeded RNG (use `rand::SeedableRng` with a fixed seed for reproducibility)
- Core gene: present in ALL k sampled genomes
- `AccumulationResult` struct implementing `DownstreamResult`

**Verification**: Unit test — known accumulation curve (e.g., 10 genomes each with unique genes) produces expected monotonic increase.

---

### Task 1.4: Add `panminer analyze` CLI subcommand

**What**: New `analyze` subcommand in `main.rs` that composes downstream analyses.

**Changes to `main.rs`**:
- Add `Analyze` variant to `Commands` enum
- Add new `analyze` match arm that:
  1. Reads `--input` (output directory from a prior `panminer` run)
  2. For each `--gwas`, `--panstripe`, `--amr`, `--neighborhood`, `--accumulation`, `--export-grapetree` flag:
     - Detect if the required tool is available
     - If available: instantiate the runner and call `run()`
     - If not available: log a warning with install instructions
  3. Write outputs to `<input>/downstream/`

**CLI arguments for `Commands::Analyze`**:
```rust
Analyze {
    /// PanMiner output directory from a prior run
    #[arg(long)]
    input: PathBuf,

    /// Run GWAS: pyseer (default), scoary2, or spydrpick
    #[arg(long)]
    gwas: bool,
    #[arg(long)]
    gwas_tool: String,  // default "pyseer"

    /// Phenotypes file (TSV: genome_id\\tphenotype)
    #[arg(long)]
    phenotypes: Option<PathBuf>,

    /// Run Panstripe evolutionary model
    #[arg(long)]
    panstripe: bool,
    /// Phylogenetic tree (Newick)
    #[arg(long)]
    tree: Option<PathBuf>,

    /// Run AMRFinderPlus resistome analysis
    #[arg(long)]
    amr: bool,
    #[arg(long)]
    amr_database: Option<PathBuf>,
    #[arg(long)]
    organism: Option<String>,

    /// Extract gene neighborhood
    #[arg(long)]
    neighborhood: bool,
    #[arg(long)]
    seed_gene: Option<String>,
    #[arg(long)]
    neighborhood_depth: Option<usize>,

    /// Generate gene accumulation curves
    #[arg(long)]
    accumulation: bool,
    #[arg(long)]
    num_samples: Option<usize>,

    /// Export for GrapeTree/iTOL visualization
    #[arg(long)]
    export_grapetree: bool,
    #[arg(long)]
    export_itol: bool,
}
```

**Verification**: `cargo build && cargo run -- analyze --help` shows the new subcommand.

---

### Task 1.5: Implement `GrapeTreeExportRunner`

**What**: Export pangenome data in GrapeTree and iTOL formats. GrapeTree invocation is optional (profiles always exported).

**File**: `src/downstream/exploration/grapetree.rs` (new file)

**Algorithm**:
1. Read `final_graph.gml` → extract genome IDs and cluster presence
2. Build allelic profile matrix: each genome = sample, each gene cluster = locus (1=present, 0=absent)
3. Write `grapetree_profiles.tsv` — tab-separated, GrapeTree input format (Newick labels + profile columns)
4. Generate `itol_annotations.txt` — iTOL dataset file:
   - Color strip dataset: core (red) vs accessory (teal) per genome
   - Binary presence heatmap for top N most variable genes
5. Optionally run GrapeTree if installed: `grapetree -i profiles.tsv -o output_prefix`
6. Write `itol_annotations.txt`

**Implementation notes**:
- GrapeTree profiles: genome name in first column, then 0/1 for each gene cluster
- iTOL format: standard iTOL dataset file with `DATASET STIPED`, `DATASET SIMPLEBAR`, etc.
- `GrapetreeResult` struct implementing `DownstreamResult`
- Tool detection: `which::which("grapetree").is_ok()`

**Verification**: Unit test — generate profiles and iTOL annotations from known graph, check TSV column counts and iTOL header format.

---

### Task 1.6: Update `lib.rs` to export `downstream` module

**What**: Add `pub mod downstream;` and re-exports to `lib.rs`.

**Changes to `lib.rs`**:
```rust
pub mod downstream;
// then downstream submodules re-exported through downstream/mod.rs
```

**Verification**: `cargo check`; `cargo test` still passes.

---

## Phase 2: GWAS + Evolutionary Tools

### Task 2.1: Relocate `PyseerRunner` to `src/downstream/gwas/pyseer.rs`

**What**: Move `src/gwas/pyseer.rs` → `src/downstream/gwas/pyseer.rs`. Keep backward-compatible re-export in `src/gwas/pyseer.rs`.

**Changes**:
1. Create `src/downstream/gwas/mod.rs` with `pub mod pyseer; pub mod scoary; pub mod spydrpick;` and re-exports
2. Create `src/downstream/gwas/pyseer.rs` (copy of existing file + enhancements)
3. Modify original `src/gwas/pyseer.rs` to just re-export from `downstream::gwas::pyseer`
4. Implement `DownstreamRunner` for `PyseerRunner`:
   - `required_inputs()` returns `[FinalGraph, PresenceAbsenceCsv, ProteinFasta, DnaFasta, PhenotypesFile]`
   - `run()` reads from output_dir, calls existing `GWASRunner::run()` internally
5. Add phenotype file builder: allow user-supplied `--phenotypes` TSV, or generate from genome metadata

**Key enhancement**: `PyseerRunner` now accepts `phenotypes_file` directly rather than generating one from `total_genes`.

**Verification**: `cargo check`; existing GWAS tests still pass.

---

### Task 2.2: Implement `Scoary2Runner`

**What**: Gene-trait association testing with Scoary2. New file `src/downstream/gwas/scoary.rs`.

**Algorithm**:
1. `is_available()`: `which::which("scoary2").is_ok()`
2. `run(output_dir)`:
   - Read `gene_presence_absence.csv` from output_dir
   - Read user-provided `phenotypes` file (or error if not provided)
   - Write temp P/A matrix in Scoary2 format (simple CSV with headers)
   - Run `scoary2 -t <phenotypes> -p <pa_matrix> -o <temp_dir>`
   - Parse `temp_dir / results.csv`: gene, trait, p_value, fdr, effect_direction, etc.
   - Write `scoary_results.csv` to `<output_dir>/downstream/`
3. `Scoary2Result::write_to()`: writes the parsed results + summary stats

**Implementation notes**:
- Scoary2 input format: `--phenotypes` TSV (genome\\tphenotype), `--presence-absence` CSV (from PanMiner)
- Output is written to a results directory: `scoary_results_<trait>.csv`
- Parse the most relevant columns: gene_id, p_value, FDR, effect_size, n_present, n_absent
- `Scoary2Result` struct implementing `DownstreamResult`

**Verification**: Unit test with synthetic P/A matrix and phenotypes — parse mock Scoary2 output.

---

### Task 2.3: Implement `SpydrPickRunner`

**What**: Gene co-selection / epistasis analysis. New file `src/downstream/gwas/spydrpick.rs`.

**Algorithm**:
1. `is_available()`: `which::which("spydrpick").is_ok()`
2. `run(output_dir)`:
   - Read `gene_presence_absence.Rtab` from output_dir (binary P/A matrix)
   - Write temp TSV in SpydrPick format
   - Run `spydrpick -i <pa_tsv> -o <temp_prefix>`
   - Parse output (SpydrPick outputs: gene pairs, MI score, p-value)
   - Write `spydrpick_correlations.csv` to `<output_dir>/downstream/`
3. `SpydrPickResult` struct implementing `DownstreamResult`

**Verification**: Unit test — parse mock SpydrPick output format.

---

### Task 2.4: Implement `PanstripeRunner`

**What**: Gene gain/loss rate estimation via GLM. New file `src/downstream/evolution/panstripe.rs`.

**Algorithm**:
1. `is_available()`: check `which::which("Rscript").is_ok()` AND R package `panstripe` is installed
2. `run(output_dir)`:
   - Read `gene_presence_absence.Rtab` from output_dir
   - Read phylogenetic tree from user-provided `--tree` flag
   - Generate an R script that calls `panstripe::run_panstripe()`:
     ```r
     library(panstripe)
     pa <- read.delim("gene_presence_absence.Rtab", row.names=1)
     tree <- read.tree("tree.nwk")
     result <- panstripe(pa, tree)
     write.csv(coef(result), "panstripe_rates.csv")
     ```
   - Write R script to temp file, run via `Rscript temp_script.R`
   - Parse `panstripe_rates.csv` output (gain_rate, loss_rate, alpha, convergence)
   - Write `panstripe_rates.csv` to `<output_dir>/downstream/`
3. `PanstripeResult` struct implementing `DownstreamResult`

**Implementation notes**:
- R package check: run `Rscript -e 'library(panstripe)'` and check exit code
- The R script can be embedded as a string in the Rust binary
- Tree file is user-provided via `--tree` flag
- `PanstripeResult::write_to()` writes `panstripe_rates.csv` and optionally copies diagnostic PDFs

**Verification**: Unit test — generate R script with known inputs, verify script syntax.

---

### Task 2.5: Implement `iTOL` phylogenetic tree annotation export

**What**: Export phylogenetic tree annotations for iTOL (as part of GrapeTreeExportRunner from Task 1.5, or as a separate lightweight runner).

**File**: `src/downstream/exploration/itol.rs` (new file, or integrated into `grapetree.rs`)

**Algorithm**:
1. Read phylogenetic tree (Newick) from user-provided `--tree` path
2. Read genome metadata from `gene_data.csv` or `final_graph.gml`
3. Generate iTOL dataset files:
   - `itol_tree.txt` — the tree in Newick format (simple pass-through)
   - `itol_annotations.txt` — color strip + text labels for each genome based on source file, cluster membership, core/accessory, etc.
4. Write to `<output_dir>/downstream/itol/`

**Implementation notes**:
- iTOL native tree format accepts Newick directly — only need to copy tree file to output
- iTOL annotation files: standard dataset format with `DATASET` header lines
- Color schemes: core=RED (#FF0000), accessory=TEAL (#008080), paralog=ORANGE (#FFA500)

**Verification**: Unit test — generate iTOL annotation file, check header format and column counts.

---

## Phase 3: Resistome + Enhancements

### Task 3.1: Implement `AmrFinderRunner`

**What**: AMR gene detection via AMRFinderPlus. New file `src/downstream/resistome/amrfinder.rs`.

**Algorithm**:
1. `is_available()`: `which::which("amrfinder").is_ok()`
2. `run(output_dir)`:
   - Read `combined_protein_CDS.fasta` and `gene_data.csv` from output_dir
   - Write protein sequences to temp input directory
   - Run `amrfinder_plus -i <input_dir> -o <output> --plus --gff <gene_data.gff>` (Plus mode = stress/virulence too)
   - Optionally add `--organism <species>` for taxon-specific point mutations
   - Parse TSV output: gene_name, scope, target_type, method, identity, coverage, contig, start, end, strand, annotation, product, resistance
   - Write `amr_results.tsv` to `<output_dir>/downstream/`
   - Generate human-readable `amr_summary.txt`: group by drug class, count genes, report evidence levels
3. `AmrFinderResult` struct implementing `DownstreamResult`

**Implementation notes**:
- GFF generation: convert `gene_data.csv` columns to GFF3 format for AMRFinderPlus
- AMRFinderPlus input: protein FASTA or assembled contigs + GFF
- Evidence levels: EXACT > ALLELE > BLAST > HMM > PARTIAL > POINT > INTERNAL_STOP
- Summary groups by `class` (drug family): beta-lactam, aminoglycoside, tetracycline, etc.
- Tool detection path: `which::which("amrfinder").or_else(|| which::which("amrfinder_plus"))`

**Verification**: Unit test — parse mock AMRFinderPlus TSV output, check evidence level parsing.

---

### Task 3.2: Scoary2 unit test with mock output

**What**: Comprehensive unit tests for all Phase 1-2 runners.

**Tests to add** (across runner files):
- `test_scoary2_parse_output` — parse known Scoary2 CSV format
- `test_spydrpick_parse_output` — parse known SpydrPick output format
- `test_panstripe_r_script_generation` — generate and validate R script syntax
- `test_amrfinder_parse_tsv` — parse known AMRFinderPlus TSV format
- `test_neighborhood_bfs_depth` — BFS finds expected nodes at each depth
- `test_accumulation_heaps_law_fit` — known data produces expected alpha
- `test_grapetree_profiles_format` — output TSV has correct column count and values
- `test_itol_annotation_format` — iTOL dataset file has correct header

---

### Task 3.3: Update `SPECS.md` and `COMPARISON.md`

**What**: Document the new downstream analysis capabilities.

**Changes**:
1. `Specs.md` — Add downstream analysis section with module structure, runner list, and CLI reference
2. `Comparison.md` — Update "What PanMiner Does Better" table with new downstream features; update Priority Roadmap to mark downstream items as complete

---

## Implementation Order

| Task | File | Type | Dependency |
|------|------|------|------------|
| **Phase 1** | | | |
| 1.1 | `src/downstream/{mod.rs,traits.rs,*/mod.rs}` | New | None |
| 1.2 | `src/downstream/exploration/neighborhood.rs` | New | 1.1 |
| 1.3 | `src/downstream/exploration/accumulation.rs` | New | 1.1 |
| 1.4 | `main.rs` (Analyze command) | Edit | 1.1 |
| 1.5 | `src/downstream/exploration/grapetree.rs` | New | 1.1 |
| 1.6 | `lib.rs` | Edit | 1.1-1.5 |
| **Phase 2** | | | |
| 2.1 | `src/downstream/gwas/pyseer.rs` (relocated) | Move | 1.1 |
| 2.2 | `src/downstream/gwas/scoary.rs` | New | 1.1 |
| 2.3 | `src/downstream/gwas/spydrpick.rs` | New | 1.1 |
| 2.4 | `src/downstream/evolution/panstripe.rs` | New | 1.1 |
| 2.5 | `src/downstream/exploration/itol.rs` | New | 1.1 |
| **Phase 3** | | | |
| 3.1 | `src/downstream/resistome/amrfinder.rs` | New | 1.1 |
| 3.2 | Unit tests across all runners | Test | 2.1-3.1 |
| 3.3 | Update SPECS.md, Comparison.md | Docs | All above |

---

## Key Design Decisions

1. **Input source**: All runners read from `<output_dir>/downstream/` files on disk (not in-memory graph/matrix). This enables the `panminer analyze <output_dir>` UX.
2. **Backward compatibility**: `src/gwas/pyseer.rs` re-exports from new location, so no API breakage.
3. **Parallel within phase**: Tasks within a phase are independent and can be done in parallel.
4. **Subprocess pattern**: All external tool wrappers follow the existing `ClipKitRunner`/`MacseRunner` pattern (detect → run → parse → write).
5. **Trait simplicity**: `DownstreamRunner::run(output_dir)` returns an opaque `Self::Output` — runners own their output data and implement `write_to()` themselves.
6. **No workflow engine**: The CLI composes runners via flags. No DAG, no dependency resolution.