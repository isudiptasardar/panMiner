# Panaroo Feature Parity — Design Spec

**Date**: 2026-04-19
**Status**: Approved — Implementation plans written
**Scope**: Address all 12 missing Panaroo features + 2 partial features

---

## Problem Statement

PanMiner is missing 12 features that Panaroo has and 2 features where our implementation is partial. These gaps affect:

- **Accuracy**: Iterative multi-threshold collapsing (#9), length-based filtering (#10), consensus removal (#5), multi-centroid nodes (#11), per-gene contig-end tracking (#13), shortest path paralog resolution (#15)
- **Completeness**: IMG/FMG parameter estimation (#1), integrate mode (#2), GFF output (#3), gene extraction (#4), standalone MSA (#7), Prodigal gene calling (#8)
- **Visualization**: Abundance plots (#6)

## Research Summary

| Gap | State-of-the-Art (2024–2025) | Chosen Approach |
|-----|------------------------------|-----------------|
| IMG/FMG | [Pangrowth](https://github.com/gi-bielefeld/pangrowth) v1.0.0 — exact pangenome growth/core curves, Heaps' law α fitting, Hill numbers. Published PCI Comp Biol 2024. | Subprocess wiring |
| Integrate mode | [PanTA](https://github.com/amromics/panta) (Genome Biology 2024) — progressive pangenome construction. 30–45× fewer CPU-hours than rebuild. | Custom implementation |
| GFF output | [PPanGGOLiN](https://ppanggolin.readthedocs.io/) `write_genomes --gff` — reconstructs GFF3 with partition/family annotations | Custom GffWriter |
| Gene extraction | Panaroo's `extract-gene` — simple FASTA lookup by cluster ID | Subcommand |
| Consensus removal | Panaroo's `remove_by_consensus` — delete spurious nodes when refound hits exceed original size | Algorithm change |
| Abundance viz | [PPanGGOLiN](https://ppanggolin.readthedocs.io/) U-shape plot, [panstripe](https://github.com/gtonkinhill/panstripe) weighted U-plot, [PangyPlot](https://pangyplot.readthedocs.io/) browser-based | D3.js HTML (matching qc_viz.rs) |
| Standalone MSA | Panaroo's `panaroo-msa` — post-run alignment command | Subcommand |
| Prodigal | [Orphos](https://github.com/FullHuman/orphos) v0.1.0 — full Rust Prodigal reimplementation, 100% compatible, Rayon-parallel, orphos-core crate on crates.io | Feature-gated crate dependency |
| Iterative collapsing | MMseqs2 cascaded clustering natively handles multi-threshold clustering with GPU support. CPU fallback loops over thresholds. | MMseqs2 cascaded + CPU iterative |
| Length filtering | Standard coverage/identity filtering in clustering (CD-HIT's `-s`/`-aL`/`-aS` flags). MMseqs2 has `--cov-mode` and `-c` for coverage. | Config parameter + both clusterers |
| Multi-centroid | Panaroo retains all centroids after merging via `centroid: list[str]`. Standard graph merging practice. | Vec<Sequence> in Node |
| Per-gene contig-end | Panaroo's `hasEnd` tracks per-gene-member contig-end status. | HashSet<GenomeId> in Node |
| Shortest path | Panaroo uses `nx.shortest_path_length` as primary paralog resolution. petgraph provides equivalent. | petgraph::algo shortest path |

---

## Architecture

### New Modules

```
src/
├── io/
│   └── orphos.rs           # Orphos/Prodigal gene calling (feature-gated: prodigal)
├── downstream/
│   └── evolution/
│       └── pangrowth.rs     # Pangrowth pangenome openness subprocess
├── output/
│   ├── gff.rs              # GFF3 output writer
│   └── abundance_viz.rs    # Gene frequency U-shape + rarefaction HTML
└── main.rs                  # New subcommands: integrate, extract-gene, msa
```

### Modified Modules

| File | Changes |
|------|---------|
| `src/config.rs` | Add `len_dif_percent`, `collapse_thresholds`, `remove_by_consensus`, `prodigal` feature |
| `src/graph/types.rs` | `Node.centroid_sequence` → `centroid_sequences: Vec<Sequence>`, `is_contig_end: bool` → `contig_end_genomes: HashSet<GenomeId>` |
| `src/graph/builder.rs` | Populate per-genome contig-end tracking |
| `src/graph/concurrent.rs` | Update merge_nodes for multi-centroid |
| `src/graph/merge.rs` | Update merge for multi-centroid |
| `src/correction/fragment.rs` | Iterative multi-threshold collapsing |
| `src/correction/missing.rs` | Add consensus removal logic |
| `src/correction/contig_end.rs` | Use per-genome contig-end data |
| `src/correction/paralog.rs` | Add shortest path resolution |
| `src/clustering/cpu.rs` | Add length-based filtering |
| `src/clustering/mmseqs.rs` | Add `--cov-mode` and `-c` flags |
| `src/pipeline.rs` | Wire iterative collapsing, consensus removal, Orphos |
| `src/output/mod.rs` | Register GffWriter, AbundanceVizWriter |
| `src/output/graph.rs` | Serialize multi-centroid and per-genome contig-end |
| `src/lib.rs` | Add new modules |
| `Cargo.toml` | Add `orphos-core` optional dep, `prodigal` feature |
| `src/main.rs` | Add Integrate, ExtractGene, Msa subcommands |

---

## Detailed Designs

### 1. Pangrowth Integration (Gap #1)

**Module**: `src/downstream/evolution/pangrowth.rs`

```rust
pub struct PangrowthRunner {
    pangrowth_path: PathBuf,
    kmer_size: usize,
    threads: usize,
}

pub struct PangrowthResult {
    pub alpha: f64,           // Heaps' law exponent
    pub kappa: f64,           // Heaps' law coefficient
    pub is_open: bool,        // alpha > 0
    pub growth_curve: Vec<(usize, usize)>,  // (n_genomes, expected_pangenome_size)
    pub core_curve: Vec<(usize, usize)>,    // (n_genomes, expected_core_size)
    pub hill_numbers: Option<HillNumbers>,
}
```

- Subprocess wrapping `pangrowth hist`, `pangrowth growth`, `pangrowth core`, `pangrowth hill`
- Input: FASTA files or presence/absence matrix
- Wired into `panminer analyze --pangrowth`
- Falls back gracefully if `pangrowth` not installed

### 2. Integrate Mode (Gap #2)

**Subcommand**: `panminer integrate`

```
panminer integrate --graph <existing_pangenome_dir> --input <new.gff> -o <output_dir>
```

**Algorithm**:
1. Load `final_graph.gml` from existing pangenome directory
2. Parse new GFF file with `GffParser`
3. Cluster new genes against existing centroids:
   - MMseqs2: Create a centroid database from the existing graph, run `easy-search` against it
   - CPU: Compute pairwise identity between new gene sequences and existing centroids
4. Add matching genes to existing nodes (increment support, add gene_members)
5. Create new nodes for genes with no match above identity threshold
6. Build adjacency edges for new genes on the same contig
7. Run correction passes on the modified graph (paralog resolution, contamination removal, contig-end pruning, fragment merging, missing gene recovery, misassembly edge cleaning)
8. Write updated output files

**Data flow**: Reuses `merge_pangenomes()` logic but for single-genome addition. The key difference is that merging combines two complete graphs, while integrating adds individual genes.

### 3. GFF Output Writer (Gap #3)

**Module**: `src/output/gff.rs`

```rust
pub struct GffWriter;

impl GffWriter {
    /// Generate GFF3 files for each genome from the corrected pangenome graph.
    /// Each gene gets its cluster_id as a Name attribute, with corrected annotations.
    pub fn write_gff(
        graph: &PangenomeGraph,
        output_dir: &Path,
    ) -> Result<Vec<PathBuf>>
}
```

**Output**: One GFF3 file per genome (`<genome_id>.gff`) with:
- `##gff-version 3` header
- `##FASTA` section with contig sequences
- CDS features with corrected cluster IDs and annotations
- `gene_data.csv` cross-references

### 4. Gene Extraction (Gap #4)

**Subcommand**: `panminer extract-gene`

```
panminer extract-gene -i <pangenome_dir> --cluster <cluster_id> -o output.fasta
```

**Implementation**:
1. Load `final_graph.gml` and `gene_data.csv`
2. Find node by cluster_id
3. Output all member sequences (DNA and/or protein) as FASTA

Simple lookup — no algorithmic complexity.

### 5. Consensus Removal in Refinding (Gap #5)

**Module**: `src/correction/missing.rs`

Add `remove_by_consensus: bool` field to `MissingGeneRecoverer`:

```rust
pub struct MissingGeneRecoverer {
    min_identity: f32,
    search_window: usize,
    prop_match: f32,
    remove_by_consensus: bool,  // NEW: delete nodes where refound hits exceed original size
}
```

**Logic**: After adding refound genes, for each node, if total refound hits for a node exceed the node's original support count, mark the node for removal. This catches spurious clusters that are artifacts of annotation errors.

**Wiring**: `CorrectionMode::Strict` enables consensus removal; `Default` and `Sensitive` disable it. Matches Panaroo's mode defaults.

### 6. Abundance Visualization (Gap #6)

**Module**: `src/output/abundance_viz.rs`

```rust
pub struct AbundanceVizWriter;

impl AbundanceVizWriter {
    /// Generate HTML report with:
    /// - Gene frequency U-shape plot (x: #genomes, y: #gene families)
    /// - Rarefaction curves (x: #genomes added, y: cumulative pangenome size)
    /// - Heaps' law fit overlay
    /// - Core/soft-core/shell/cloud partition bars
    pub fn write_abundance_report(
        matrix: &BitPackedMatrix,
        heaps_fit: &HeapsLawFit,
        output_path: &Path,
    ) -> Result<()>
}
```

**Approach**: HTML with embedded D3.js (matching `qc_viz.rs` pattern). No Python/matplotlib dependency. Wired into `--formats abundance` or `panminer analyze --abundance`.

### 7. Standalone MSA (Gap #7)

**Subcommand**: `panminer msa`

```
panminer msa -i <pangenome_dir> --aligner mafft --mode core -o alignments/
```

**Implementation**: Load `final_graph.gml` + `gene_data.csv`, identify core/accessory genes, run MSA using existing `AlignmentRunner` trait. Reuses `MafftRunner`, `ClustalOmegaRunner`, `PrankRunner`.

### 8. Orphos/Prodigal Gene Calling (Gap #8)

**Module**: `src/io/orphos.rs` (feature-gated: `prodigal`)

```rust
pub struct OrphosRunner {
    config: OrphosConfig,
}

impl OrphosRunner {
    /// Run Prodigal/Orphos gene prediction on a FASTA file.
    /// Returns a list of predicted genes with coordinates.
    pub fn predict_genes(&self, fasta_path: &Path) -> Result<Vec<PredictedGene>>

    /// Run in metagenomic mode.
    pub fn predict_genes_metagenomic(&self, fasta_path: &Path) -> Result<Vec<PredictedGene>>
}
```

**Dependency**: `orphos-core = { version = "0.1.0", optional = true }` behind `prodigal` feature flag. No subprocess needed — pure Rust library call.

**Pipeline wiring**: When `--pipeline-mode prodigal` is set, use `OrphosRunner` for gene calling on unannotated FASTA inputs instead of requiring pre-annotated GFF3.

### 9. Iterative Multi-Threshold Collapsing (Gap #9)

**Module**: `src/correction/fragment.rs` (modified)

**Current behavior**: Single-pass at `collapse_threshold: 0.70`.

**New behavior**:

```rust
pub struct FragmentMerger {
    // ...existing fields...
    collapse_thresholds: Vec<f32>,  // NEW: default [0.99, 0.95, 0.9, 0.8, 0.7]
}
```

**MMseqs2 path**: Use MMseqs2 cascaded clustering with `--cluster-reuse` flag, which naturally handles progressive clustering at decreasing thresholds in a single GPU-accelerated pass. The thresholds are passed via `--cluster-steps` parameter.

**CPU fallback path**: Loop over `collapse_thresholds`, calling `collapse_gene_families_with_cache` at each step. The `DistanceCache` is preserved across iterations, matching Panaroo's reuse of the distance matrix.

**Threshold order**: Iterate from high to low identity (e.g., [0.99, 0.95, 0.9, 0.8, 0.7]). At each threshold, only merge pairs not already merged at a higher threshold. This matches Panaroo's behavior where more stringent merges happen first.

**Pipeline wiring**: Replace the two current `collapse_gene_families_with_cache` calls with a loop over `collapse_thresholds`. After mistranslation correction (which stays at identity 0.99), run the family collapse loop.

### 10. Length-Based Filtering During Clustering (Gap #10)

**Module**: `src/clustering/cpu.rs` and `src/clustering/mmseqs.rs`

**Config**: Add to `PanminerConfig`:
```rust
/// Length difference cutoff for clustering (0.0–1.0, default 0.98).
/// Gene pairs with length difference > (1 - len_dif_percent) are excluded.
pub len_dif_percent: f32,
```

**CPU clusterer**: In `CpuClusterer::cluster()`, after computing sequence identity, check:
```rust
let max_len = a.len().max(b.len()) as f32;
let len_diff = (a.len().abs_diff(b.len())) as f32 / max_len;
if len_diff > (1.0 - self.len_dif_percent) {
    continue; // skip this pair
}
```

**MMseqs2**: Pass `--cov-mode 1 -c {len_dif_percent}` to `easy-cluster`.

### 11. Multi-Centroid Nodes (Gap #11)

**Module**: `src/graph/types.rs`

**Change**: `Node.centroid_sequence: Option<Sequence>` → `centroid_sequences: Vec<Sequence>`

```rust
pub struct Node {
    // ...existing fields...
    /// Centroid sequences (one per original cluster; multiple after merging)
    pub centroid_sequences: Vec<Sequence>,
    // ...existing fields...
}
```

**Impact sites**:
- `GeneCluster.centroid: Option<Sequence>` → `centroids: Vec<Sequence>`
- `GraphBuilder`: populate from cluster centroids
- `ConcurrentGraph::merge_nodes`: concatenate centroids instead of picking one
- `FragmentMerger`: use all centroids for identity comparison
- `MissingGeneRecoverer`: compare against all centroids
- `GmlWriter`: serialize all centroids
- `GmlWriter` reading: parse multi-centroid GML
- `OutputWriter`: use first centroid as representative for FASTA output, all centroids for reference
- `JsonWriter`: serialize all centroid sequences
- `MatrixWriter`: use first centroid for representative sequence

**Migration**: A `Node` with a single centroid becomes `centroid_sequences: vec![seq]`. Empty becomes `centroid_sequences: vec![]`. The `Option<Sequence>` → `Vec<Sequence>` change is non-breaking at the data level since `vec![]` is the `None` equivalent.

**GML backward compatibility**: When reading GML files, check for the old `centroid_sequence` attribute. If found, wrap into `vec![seq]`. If the new `centroid_sequences` attribute is found (JSON array), parse normally. When writing GML, always use the new `centroid_sequences` attribute in JSON format.

### 12→13. Per-Gene Contig-End Tracking (Gaps #12, #13)

**Module**: `src/graph/types.rs`

**Change**: `Node.is_contig_end: bool` → `contig_end_genomes: HashSet<GenomeId>`

```rust
pub struct Node {
    // ...existing fields...
    /// Genomes where this cluster appears at a contig end.
    /// Replaces the boolean is_contig_end for per-gene tracking.
    pub contig_end_genomes: HashSet<GenomeId>,
    // ...existing fields...
}
```

**Impact sites**:
- `GraphBuilder`: When marking contig-end genes, add the genome_id to `contig_end_genomes` instead of setting a bool
- `ContigEndPruner`: Change condition from `node.is_contig_end` to `!node.contig_end_genomes.is_empty()`. For removal, only remove if ALL member genomes are in `contig_end_genomes` (or use proportion-based threshold)
- `MisassemblyEdgeCleaner`: Change `hasEnd` check to `!node.contig_end_genomes.is_empty()`
- `GmlWriter`: Serialize as comma-separated genome IDs
- `GmlWriter` reading: Parse back into HashSet

**Backward compatibility**: When reading old GML files that have `is_contig_end true/false`, convert `true` → all member genomes, `false` → empty set.

### 15. Shortest Path Paralog Resolution (Gap #15)

**Module**: `src/correction/paralog.rs`

**Change**: Add shortest path as primary resolution method, context vectors as fallback.

```rust
impl ParalogResolver {
    fn resolve_paralogs(&self, graph: &mut ConcurrentGraph) -> Result<ParalogStats> {
        // For each paralog group:
        // 1. Try shortest path distance between paralog copies
        // 2. If no path exists within max_context depth, fall back to context vector similarity
        // 3. Merge paralog copies into the best-matching cluster
    }

    fn shortest_path_distance(
        graph: &ConcurrentGraph,
        from: &ClusterId,
        to: &ClusterId,
    ) -> Option<usize> {
        // Use petgraph::algo::bfs to find shortest path
        // Returns None if no path within max_context depth
    }
}
```

**Implementation**: Use `petgraph::algo::bfs` or custom BFS since `ConcurrentGraph` uses `DashMap` not `petgraph::Graph`. The BFS is simple: start from `from`, explore neighbors up to `max_context` depth, return distance if `to` is reached.

---

## Configuration Changes

### `PanminerConfig` Additions

```rust
// Clustering
pub len_dif_percent: f32,                    // default: 0.98
pub collapse_thresholds: Vec<f32>,           // default: [0.99, 0.95, 0.9, 0.8, 0.7]

// Correction
pub remove_by_consensus: Option<bool>,        // None = mode-dependent

// Pipeline
pub prodigal_mode: bool,                     // --pipeline-mode prodigal
```

### CLI Additions

```
# Existing command additions
--len-dif-percent 0.98          # Length difference cutoff
--collapse-thresholds 0.99,0.95,0.9,0.8,0.7  # Iterative collapsing thresholds
--remove-by-consensus           # Enable consensus removal (strict default)

# New subcommands
panminer integrate --graph <dir> --input <new.gff> -o <dir>
panminer extract-gene -i <dir> --cluster <id> -o output.fasta
panminer msa -i <dir> --aligner mafft --mode core -o alignments/
panminer analyze --pangrowth    # Pangrowth pangenome openness
```

### Feature Flag Addition

```toml
[features]
prodigal = ["orphos-core"]
```

---

## Implementation Order

The implementation is ordered by dependency chain and impact:

| Phase | Gaps | Rationale |
|-------|------|-----------|
| **P0** | #11, #12→13 | Data structure changes that other features depend on |
| **P1** | #9, #10, #5, #15 | Core algorithm changes (iterative collapsing, length filtering, consensus removal, shortest path) |
| **P2** | #2, #3, #4, #7 | New subcommands (integrate, GFF output, extract-gene, MSA) |
| **P3** | #1, #6, #8 | External tool wiring (Pangrowth, abundance viz, Orphos) |

Within each phase, features are independent and can be implemented in parallel.

---

## Testing Strategy

Each gap gets:
- **Unit tests**: Testing the new behavior in isolation
- **Integration tests**: Testing the full pipeline with the new features
- **Regression tests**: Ensuring existing tests still pass

Key test scenarios:
- P0: Multi-centroid merge preserves all sequences; per-genome contig-end tracking correct
- P1: Iterative collapsing produces more merged clusters than single-pass; length filtering rejects pairs with >2% length difference; consensus removal deletes spurious nodes; shortest path resolves paralogs that context vectors miss
- P2: Integrate adds a genome to existing graph; GFF output is valid GFF3; extract-gene produces correct FASTA; standalone MSA works
- P3: Pangrowth produces openness classification; abundance HTML renders correctly; Orphos predicts genes from FASTA

---

## Risk Assessment

| Risk | Mitigation |
|------|-----------|
| Multi-centroid change touches many files | Phase P0 first; comprehensive grep for all `centroid_sequence` references |
| Iterative collapsing may change results significantly | Make `collapse_thresholds` configurable; default matches Panaroo |
| Orphos crate is new (v0.1.0, GPL-3.0) | Feature-gated; not required for core functionality; GPL compatible with our project |
| Pangrowth binary dependency | Graceful fallback with clear error message if not installed |
| Per-genome contig-end changes correction behavior | Run existing test suite; compare results with boolean approach |
| Integrate mode is complex | Start simple (add genes, run corrections); defer incremental correction optimization |