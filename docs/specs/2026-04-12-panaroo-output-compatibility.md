# Panaroo Output Compatibility Design

**Date:** 2026-04-12
**Status:** Approved
**Depends on:** cDBG backbone (docs/specs/2026-04-12-cdbg-backbone-design.md) — independent but complementary

## Context

PanMiner produces output files that are structurally similar to Panaroo but have 5 compatibility gaps that prevent downstream tools (pyseer, scoary, Panaroo visualizers) from working unmodified. The root cause is that PanMiner's `Node` struct discards per-gene detail when converting `GeneCluster` to `Node`, blocking gene-member-aware output formats.

## Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Data model | Enrich `Node` with `gene_members: HashMap<GenomeId, Vec<String>>` | Unblocks all 5 gaps; single root-cause fix |
| Roary CSV | Always write alongside `gene_presence_absence.csv` | Matches Panaroo behavior; tools expect `_roary.csv` filename |
| gene_data.csv | Add `dna_sequence` + `protein_sequence` columns | Panaroo includes sequences; data already available from `centroid_sequence` |
| GML attributes | Expand to include `length`, `seq`, `protein`, `genome_ids`, `member` | Panaroo GML is parsed by downstream visualizers |
| BMGE filtering | Python subprocess via Biopython | BMGE is a Python library, not a standalone binary; matches our subprocess pattern |
| pre_filt_graph.gml | Track in `OutputPaths` | Already written in pipeline; just needs formal tracking |

## Gap 1: `gene_presence_absence_roary.csv`

### Problem
Panaroo produces two P/A files:
- `gene_presence_absence.csv` (3 metadata columns)
- `gene_presence_absence_roary.csv` (14 metadata columns, semicolon-delimited gene IDs per genome)

PanMiner writes one file with 14 columns but uses the simplified filename. Downstream tools (scoary, pyseer) look for `_roary.csv`.

### Fix
- Add `write_roary_compat_csv()` method to `MatrixWriter` in `src/output/matrix.rs`
- Per-genome cells: semicolon-delimited gene IDs from `node.gene_members[genome_id]`
- e.g., `geneA;geneB` instead of just `cluster_42`
- Also write the simplified 3-column `gene_presence_absence.csv` (Gene, Annotation, Genome columns)
- Add `matrix_roary_csv: Option<PathBuf>` to `OutputPaths`
- Always written when `OutputFormat::Matrix` is selected

### Output Format
```
Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome Fragment,Order within Fragment,Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,Avg group size nuc,genome1,genome2,...
cluster_0,cluster_0,hypothetical protein,3,3,1.00,,,,,,450,450,450.0,geneA;geneB,,geneC
```

## Gap 2: `gene_data.csv` with DNA/protein sequences

### Problem
PanMiner writes 8 columns: `gene_id,gene_name,annotation,contig,start,end,strand,support`. Location columns are all `NA`. Panaroo includes DNA and protein sequences.

### Fix
- Add `dna_sequence` and `protein_sequence` columns after `support`
- DNA from `node.centroid_sequence`
- Protein from `translate(node.centroid_sequence)`
- Also populate `contig`, `start`, `end`, `strand` from the first gene member in `node.gene_members` (look up in the original gene data)
- New header: `gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence`
- Modify `write_gene_data()` in `src/output/json.rs`

### Implementation Note
To populate contig/start/end/strand, the graph builder needs to preserve a `HashMap<GeneId, Gene>` lookup table that output writers can access. Add this as a field on `PangenomeGraph`:
```rust
pub gene_lookup: HashMap<GeneId, Gene>,
```

## Gap 3: GML attribute expansion

### Problem
PanMiner GML nodes have only `id`, `label`, `support`, `is_paralog`, `annotation`. Panaroo includes `length`, `seq`, `protein`, `genome_ids`, `member`.

### Fix
Expand `GmlWriter::write()` in `src/output/graph.rs`:

**Node attributes to add:**
- `length` — `node.centroid_sequence.as_ref().map(|s| s.len()).unwrap_or(0)`
- `seq` — centroid DNA sequence (GML-escaped string)
- `protein` — translated centroid sequence
- `genome_ids` — comma-separated genome IDs from `node.genomes`
- `member` — semicolon-separated gene IDs from `node.gene_members.values().flatten()`

**Edge attributes to add:**
- `genome_ids` — comma-separated genome IDs from `edge.genomes`

**GML string escaping:** Replace `"` with `\"` and `\` with `\\` in string values.

## Gap 4: BMGE-filtered alignment

### Problem
Panaroo offers BMGE (Block Mapping and Gathering with Entropy) for entropy-based column filtering. PanMiner only has ClipKIT trimming.

### Fix
- Add `BmgeRunner` in `src/output/trim.rs` alongside `ClipKitRunner`
- BMGE runs as a Python subprocess using Biopython:
```python
from Bio.Align import AlignInfo
from Bio.Align.Applications import BMGECommandline
```
Actually, BMGE is typically invoked via Java (`java -jar BMGE.jar`) or Biopython. For simplicity, use the Biopython approach since it's the most common install method:
```bash
python3 -c "
import sys
from Bio import AlignIO
from bmge import bmge
alignment = AlignIO.read(sys.argv[1], 'fasta')
filtered = bmge(alignment, gap_threshold=0.2)
AlignIO.write(filtered, sys.argv[2], 'fasta')
" input.aln output.aln
```
- Detect via `python3 -c "import bmge"` or check for `bmge.jar` on PATH
- Add `--filter-alignment bmge` CLI flag (alongside existing `--trim-alignment` for ClipKIT)
- Output filename: `core_gene_alignment.BMGE.aln`
- Add `bmge_alignment: Option<PathBuf>` to `OutputPaths`

## Gap 5: `pre_filt_graph.gml` tracking

### Problem
`pre_filt_graph.gml` is already written in `pipeline.rs:170` but not tracked in `OutputPaths`.

### Fix
- Add `pre_filt_graph: Option<PathBuf>` to `OutputPaths`
- Return the path from the pipeline write step
- Ensure it uses the same enriched GML attributes (once Gap 3 is resolved)

## Root Cause: Node Data Model Enrichment

### Current `Node` struct
```rust
pub struct Node {
    pub cluster_id: ClusterId,
    pub support: usize,
    pub genomes: HashSet<GenomeId>,
    pub annotations: HashSet<String>,
    pub is_paralog: bool,
    pub centroid_sequence: Option<Sequence>,
    pub is_contig_end: bool,
    pub contig_sequences: HashMap<String, Sequence>,
}
```

### Proposed additions
```rust
pub struct Node {
    // ... existing fields ...

    /// Gene members per genome: genome_id -> [gene_id, gene_id, ...]
    pub gene_members: HashMap<GenomeId, Vec<String>>,
}
```

### Changes to `Node::from_cluster()`
Current: Takes only `&GeneCluster`, which has `genes: Vec<GeneId>` but no genome mapping.

New signature:
```rust
pub fn from_cluster_with_genes(
    cluster: &GeneCluster,
    gene_data: &HashMap<GeneId, Gene>,
) -> Self
```

This populates `gene_members` by looking up each gene's `genome_id` from the `Gene` struct.

### Changes to `PangenomeGraph`
Add a gene lookup table:
```rust
pub gene_lookup: HashMap<GeneId, Gene>,
```

This is populated during graph building and enables output writers to look up contig/start/end/strand for any gene ID.

## File Changes Summary

| File | Action | Description |
|---|---|---|
| `src/graph/types.rs` | **Modify** | Add `gene_members` field to `Node`, add `from_cluster_with_genes()`, add `gene_lookup` to `PangenomeGraph` |
| `src/graph/builder.rs` | **Modify** | Pass gene data to `from_cluster_with_genes()`, populate `gene_lookup` |
| `src/output/matrix.rs` | **Modify** | Add `write_roary_compat_csv()` with semicolon-delimited gene IDs |
| `src/output/json.rs` | **Modify** | Add `dna_sequence`, `protein_sequence` columns to `gene_data.csv` |
| `src/output/graph.rs` | **Modify** | Add `length`, `seq`, `protein`, `genome_ids`, `member` node attributes; `genome_ids` edge attribute |
| `src/output/trim.rs` | **Modify** | Add `BmgeRunner` for Python/Biopython BMGE |
| `src/output/mod.rs` | **Modify** | Add `matrix_roary_csv`, `bmge_alignment`, `pre_filt_graph` to `OutputPaths` |
| `src/main.rs` | **Modify** | Add `--filter-alignment` CLI flag |
| `src/config.rs` | **Modify** | Add `filter_method` config option |
| `src/pipeline.rs` | **Modify** | Track `pre_filt_graph.gml` in output paths |

## Verification Checklist

- [ ] `gene_presence_absence_roary.csv` is written with semicolon-delimited gene IDs
- [ ] `gene_data.csv` includes `dna_sequence` and `protein_sequence` columns
- [ ] GML nodes include `length`, `seq`, `protein`, `genome_ids`, `member`
- [ ] GML edges include `genome_ids`
- [ ] `--filter-alignment bmge` produces `core_gene_alignment.BMGE.aln`
- [ ] `pre_filt_graph.gml` is tracked in `OutputPaths`
- [ ] `scoary --pa gene_presence_absence_roary.csv` works unmodified
- [ ] `pyseer --pheno phenotypes.txt gene_presence_absence_roary.csv` works