---
name: missing_gene_recovery_plan
description: Missing gene recovery implementation plan
type: project
---

**Panaroo's Approach (from docs/gettingstarted/params.md):**
- `find_missing()`: Re-finds genes missed by Prokka
- Parameters: `search_radius` (default 5000bp), `prop_match` (default 0.2)
- Strategy: Search flanking regions of missing gene locations

**PanMiner Implementation Plan:**

**Phase 1: Data Structure Changes**
1. Store contig sequences in `Node` or `GeneCluster`
2. Add `flanking_sequences: HashMap<GenomeId, Vec<u8>>` to nodes

**Phase 2: Missing Gene Detection**
1. Identify clusters present in some but not all genomes
2. For each "missing" cluster in a genome, search flanking regions

**Phase 3: K-mer Search Algorithm**
1. Extract 11-mers from cluster centroid sequence
2. For each genome missing the cluster:
   - Extract 5000bp flanking regions from contig
   - Search for k-mer matches with >=70% identity
   - If found, add gene to cluster with partial support

**Phase 4: Integration**
1. Wire into `run_corrections()` pipeline
2. Run after contamination removal, before matrix construction

**Key Parameters:**
- K-mer size: 11 (Panaroo default)
- Search radius: 5000bp
- Match threshold: 70% identity

**Web Sources:**
- Panaroo params: https://github.com/gtonkinhill/panaroo/blob/master/docs/gettingstarted/params.md
- Panaroo __main__.py: https://github.com/gtonkinhill/panaroo/blob/master/panaroo/__main__.py
