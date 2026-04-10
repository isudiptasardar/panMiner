---
name: panaroo_algorithm_reference
description: Panaroo algorithm reference for PanMiner development
type: project
---

Panaroo is a Python-based pangenome analysis tool with the following key components:

**Clustering**: Uses CD-HIT v4.8.1 at 98% identity threshold with length difference cutoff of 0.98. Protein families are clustered at a lower threshold of 0.7.

**Graph Construction**: NetworkX-based graph where nodes are gene clusters (COGs) and edges connect clusters when genes are adjacent on any contig. Paralogs are split into separate clusters by default.

**Error Correction Pipeline** (in order):
1. `trim_low_support_trailing_ends()` - Recursively removes low-support degree-1 nodes at contig ends
2. `collapse_families()` with `correct_mistranslations=True` - Merges genes with ≥95% coverage and ≥99% identity (corrects different reading frames)
3. `collapse_families()` with `correct_families=True` - Merges gene families sharing neighbors at ≥70% identity
4. `find_missing()` - K-mer search (11-mers, 5000bp flanking) to refind missing genes
5. `clean_misassembly_edges()` - Removes low-support edges flagged as misassemblies

**Correction Modes**:
- `strict`: ≥5% genome presence threshold, aggressive removal
- `moderate`: ≥1% presence, keeps refound genes
- `sensitive`: No deletion, only merge/refind operations

**Output Files** (Roary-compatible naming):
- `gene_presence_absence.csv` - Presence/absence matrix
- `gene_presence_absence.Rtab` - Binary format
- `final_graph.gml` - NetworkX graph
- `struct_presence_absence.csv` - Structural variants
- `pan_genome_reference.fa` - Reference genome
- `core_gene_alignment.aln` - Core alignment

**Input Requirements**: GFF3 (Prokka format preferred), GenBank (.gbk/.gb/.gbff). Requires Prokka re-annotation for best results.

**Key Python Libraries**: NetworkX, Biopython, edlib (for pwdist_edlib pairwise distances), subprocess calls to cd-hit.

**Thread Safety Issue**: Known KeyError in `collapse_families` with high thread counts (60 threads); workaround is 8 threads. Caused by `pwdist_edlib` sub-function under concurrent load.

**Reference**: Tonkin-Hill et al., Genome Biology (2020) - https://doi.org/10.1186/s13059-020-02090-4

**How to apply**: When implementing similar functionality in PanMiner:
- Use MMseqs2 instead of CD-HIT (GPU-capable)
- Use DashMap for concurrent graph instead of NetworkX
- Implement the same correction logic but with Rust data structures
- PanMiner currently has the correction modules but needs real sequence data passed to FragmentMerger
