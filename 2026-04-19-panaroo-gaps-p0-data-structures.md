# Panaroo Feature Parity — Phase 0: Data Structure Changes

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Change `Node.centroid_sequence: Option<Sequence>` → `centroid_sequences: Vec<Sequence>` and `Node.is_contig_end: bool` → `contig_end_genomes: HashSet<GenomeId>` to support multi-centroid nodes and per-gene contig-end tracking.

**Architecture:** These are foundational data structure changes that all other phases depend on. `centroid_sequences: Vec<Sequence>` preserves sequence diversity after node merges. `contig_end_genomes: HashSet<GenomeId>` tracks which genomes have this cluster at a contig end, enabling more precise correction.

**Tech Stack:** Rust, thiserror, serde, DashMap, HashSet

---

### Task 1: Change `GeneCluster.centroid` to `centroids: Vec<Sequence>`

**Files:**
- Modify: `src/graph/types.rs:141-161` (GeneCluster struct)
- Modify: `src/clustering/cpu.rs:74` (CpuClusterer initialization)
- Modify: `src/clustering/mmseqs.rs:193` (MMseqsRunner initialization)
- Modify: `src/graph/merge.rs:326` (merge test helper)
- Modify: `tests/fragment_merge_test.rs:14,20,63,80,106,111` (test data)
- Modify: `src/correction/paralog.rs:502,507` (paralog test data)

- [ ] **Step 1: Write failing tests for `centroids: Vec<Sequence>`**

In `src/graph/types.rs`, add a test that creates a `GeneCluster` with multiple centroids:

```rust
#[test]
fn test_gene_cluster_multiple_centroids() {
    let cluster = GeneCluster {
        id: ClusterId::new("cluster_0"),
        genes: vec![GeneId::new("gene_0")],
        centroids: vec![b"ATCG".to_vec(), b"GCTA".to_vec()],
        is_paralog: false,
        support: 1,
    };
    assert_eq!(cluster.centroids.len(), 2);
    assert_eq!(cluster.centroids[0], b"ATCG".to_vec());
    assert_eq!(cluster.centroids[1], b"GCTA".to_vec());
}
```

- [ ] **Step 2: Change `GeneCluster.centroid: Option<Sequence>` to `centroids: Vec<Sequence>` in types.rs**

In `src/graph/types.rs`, change the `GeneCluster` struct:

```rust
// Before:
pub centroid: Option<Sequence>,

// After:
pub centroids: Vec<Sequence>,
```

Change the `new()` method default:

```rust
// Before:
centroid: None,

// After:
centroids: vec![],
```

- [ ] **Step 3: Update CpuClusterer in cpu.rs**

Change `src/clustering/cpu.rs:74`:

```rust
// Before:
new_cluster.centroid = Some(gene.sequence.clone());

// After:
new_cluster.centroids = vec![gene.sequence.clone()];
```

- [ ] **Step 4: Update MMseqsRunner in mmseqs.rs**

Change `src/clustering/mmseqs.rs:193`:

```rust
// Before:
cluster.centroid = Some(seq.clone());

// After:
cluster.centroids = vec![seq.clone()];
```

- [ ] **Step 5: Update merge test helper in merge.rs**

Change `src/graph/merge.rs:326` and any test helpers that set `cluster.centroid`:

```rust
// Before:
cluster.centroid = centroid;  // where centroid: Option<Sequence>

// After:
cluster.centroids = centroid;  // where centroid: Vec<Sequence>
```

Update the `make_test_graph` function signature if it takes `Option<Sequence>` to take `Vec<Sequence>`.

- [ ] **Step 6: Update paralog test data in paralog.rs**

Change `src/correction/paralog.rs:502,507`:

```rust
// Before:
para1.centroid = Some(b"ATCGATCGATCGATCG".to_vec());

// After:
para1.centroids = vec![b"ATCGATCGATCGATCG".to_vec()];
```

Same pattern for `para2`.

- [ ] **Step 7: Update fragment_merge_test.rs**

Change all 6 occurrences of `cluster.centroid = Some(...)` to `cluster.centroids = vec![...]` and `cluster.centroid = None` to `cluster.centroids = vec![]`:

```rust
// Before:
cluster_a.centroid = Some(b"ATCGATCGATCGATCGATCGATCG".to_vec());

// After:
cluster_a.centroids = vec![b"ATCGATCGATCGATCGATCGATCG".to_vec()];
```

- [ ] **Step 8: Run all tests to verify GeneCluster changes**

Run: `cargo test --features full 2>&1 | head -80`
Expected: All tests pass. The `GeneCluster.centroids` field is now `Vec<Sequence>`.

- [ ] **Step 9: Commit GeneCluster changes**

```bash
git add src/graph/types.rs src/clustering/cpu.rs src/clustering/mmseqs.rs src/graph/merge.rs src/correction/paralog.rs tests/fragment_merge_test.rs
git commit -m "refactor: change GeneCluster.centroid to centroids Vec<Sequence>

Multi-centroid support: merged nodes retain all centroid sequences
instead of picking one. Vec<Sequence> preserves sequence diversity."
```

---

### Task 2: Change `Node.centroid_sequence` to `centroid_sequences: Vec<Sequence>`

**Files:**
- Modify: `src/graph/types.rs:207,226` (Node struct)
- Modify: `src/graph/concurrent.rs:246-253` (merge_nodes)
- Modify: `src/pipeline.rs:731,739,763,796` (pipeline references)
- Modify: `src/correction/paralog.rs:116,117` (paralog resolution)
- Modify: `src/clustering/alignment_traits.rs:81-82` (alignment builder)
- Modify: `src/graph/merge.rs:261,273` (merge centroids)
- Modify: `src/output/json.rs:117,120,147,166,185,280` (JSON output)
- Modify: `src/output/alignment.rs:50,91` (alignment output)
- Modify: `src/output/graph.rs:38,42,48,148` (GML output)
- Modify: `tests/fragment_merge_test.rs:43,67,84,94,126` (test references)

- [ ] **Step 1: Change `Node.centroid_sequence: Option<Sequence>` to `centroid_sequences: Vec<Sequence>`**

In `src/graph/types.rs`, change the `Node` struct:

```rust
// Before:
pub centroid_sequence: Option<Sequence>,

// After:
pub centroid_sequences: Vec<Sequence>,
```

Change the `from_cluster` initialization:

```rust
// Before:
centroid_sequence: cluster.centroid.clone(),

// After:
centroid_sequences: cluster.centroids.clone(),
```

Change the `from_cluster_with_genes` method similarly.

- [ ] **Step 2: Update `merge_nodes` in concurrent.rs**

In `src/graph/concurrent.rs`, add centroid merging to `merge_nodes`:

```rust
// After the existing merge lines (support, annotations, etc.), add:
target_node.centroid_sequences.extend(source_node.centroid_sequences);
target_node.contig_end_genomes.extend(source_node.contig_end_genomes);
```

- [ ] **Step 3: Update pipeline.rs references (4 locations)**

At lines 731, 739, 763, and 796, change `node.centroid_sequence.clone().unwrap_or_default()` to use the first centroid or iterate:

```rust
// For representative sequence extraction (pipeline needs one sequence per cluster):
// Before:
(cluster_id.to_string(), node.centroid_sequence.clone().unwrap_or_default())

// After:
(cluster_id.to_string(), node.centroid_sequences.first().cloned().unwrap_or_default())
```

For line 796 (cluster sequence extraction for missing gene recovery):

```rust
// Before:
node.centroid_sequence.clone().map(|seq| (cluster_id, seq))

// After:
node.centroid_sequences.first().map(|seq| (cluster_id.clone(), seq.clone()))
```

- [ ] **Step 4: Update paralog.rs references**

Change `src/correction/paralog.rs:116,117`:

```rust
// Before:
let seq_a = node_a.centroid_sequence.as_deref().unwrap_or(&[]);
let seq_b = node_b.centroid_sequence.as_deref().unwrap_or(&[]);

// After:
let seq_a = node_a.centroid_sequences.first().map(|s| s.as_slice()).unwrap_or(&[]);
let seq_b = node_b.centroid_sequences.first().map(|s| s.as_slice()).unwrap_or(&[]);
```

- [ ] **Step 5: Update alignment_traits.rs**

Change `src/clustering/alignment_traits.rs:81-82`:

```rust
// Before:
if let Some(centroid) = &node.centroid_sequence {
    sequences.push((cluster_id.to_string(), centroid.clone()));
}

// After:
if let Some(centroid) = node.centroid_sequences.first() {
    sequences.push((cluster_id.to_string(), centroid.clone()));
}
```

- [ ] **Step 6: Update merge.rs**

Change `src/graph/merge.rs:261`:

```rust
// Before:
all_centroids.push((cluster_id.clone(), node.centroid_sequence.clone()));

// After:
for seq in &node.centroid_sequences {
    all_centroids.push((cluster_id.clone(), Some(seq.clone())));
}
```

- [ ] **Step 7: Update json.rs (6 locations)**

For each `node.centroid_sequence` reference in `src/output/json.rs`, change to `node.centroid_sequences.first()`:

```rust
// Before (line 117):
let dna_seq = node.centroid_sequence.as_ref()

// After:
let dna_seq = node.centroid_sequences.first()

// Before (line 120):
let protein_seq = if let Some(seq) = &node.centroid_sequence {

// After:
let protein_seq = if let Some(seq) = node.centroid_sequences.first() {
```

Same pattern for lines 147, 166, 185.

For test at line 280:

```rust
// Before:
node.centroid_sequence = Some(b"ATGCGT".to_vec());

// After:
node.centroid_sequences = vec![b"ATGCGT".to_vec()];
```

- [ ] **Step 8: Update alignment.rs (2 locations)**

```rust
// Before:
node.centroid_sequence.clone()

// After:
node.centroid_sequences.first().cloned()
```

- [ ] **Step 9: Update graph.rs (GML output) (4 locations)**

For representative sequence in GML:

```rust
// Before (line 38):
let length = node.centroid_sequence.as_ref().map(|s| s.len()).unwrap_or(0);

// After:
let length = node.centroid_sequences.first().map(|s| s.len()).unwrap_or(0);
```

For GML serialization, write all centroids as JSON array:

```rust
// Before (lines 42, 48):
if let Some(seq) = &node.centroid_sequence {
    // write seq
}

// After:
if !node.centroid_sequences.is_empty() {
    let centroid_seqs: Vec<String> = node.centroid_sequences.iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();
    // Write as JSON array
}
```

For test at line 148:

```rust
// Before:
node.centroid_sequence = Some(b"ATGCGT".to_vec());

// After:
node.centroid_sequences = vec![b"ATGCGT".to_vec()];
```

- [ ] **Step 10: Update fragment_merge_test.rs (5 locations)**

```rust
// Before:
entry.value().centroid_sequence.clone().map(|seq| (entry.key().to_string(), seq))

// After:
entry.value().centroid_sequences.first().cloned().map(|seq| (entry.key().to_string(), seq))
```

For `node.centroid_sequence = None`:

```rust
// Before:
node.centroid_sequence = None;

// After:
node.centroid_sequences = vec![];
```

For assertions:

```rust
// Before:
assert_eq!(node.centroid_sequence, Some(b"ATCGATCGATCGATCG".to_vec()));

// After:
assert_eq!(node.centroid_sequences, vec![b"ATCGATCGATCGATCG".to_vec()]);
```

- [ ] **Step 11: Run all tests**

Run: `cargo test --features full 2>&1 | head -80`
Expected: All tests pass with `centroid_sequences: Vec<Sequence>`.

- [ ] **Step 12: Commit Node.centroid_sequence changes**

```bash
git add src/graph/types.rs src/graph/concurrent.rs src/pipeline.rs src/correction/paralog.rs src/clustering/alignment_traits.rs src/graph/merge.rs src/output/json.rs src/output/alignment.rs src/output/graph.rs tests/fragment_merge_test.rs
git commit -m "refactor: change Node.centroid_sequence to centroid_sequences Vec<Sequence>

Multi-centroid support: merged nodes now retain all centroid sequences.
GML, JSON, and alignment output use first centroid as representative."
```

---

### Task 3: Change `Node.is_contig_end: bool` to `contig_end_genomes: HashSet<GenomeId>`

**Files:**
- Modify: `src/graph/types.rs:209,227` (Node struct)
- Modify: `src/graph/builder.rs:133-138,336-338` (GraphBuilder construction + test)
- Modify: `src/graph/concurrent.rs:246-253` (merge_nodes — already updated in Task 2)
- Modify: `src/correction/contig_end.rs:49,114,150,187` (ContigEndPruner + tests)
- Modify: `src/correction/misassembly.rs:67,137` (MisassemblyEdgeCleaner + test)

- [ ] **Step 1: Change `Node.is_contig_end: bool` to `contig_end_genomes: HashSet<GenomeId>`**

In `src/graph/types.rs`, change the `Node` struct:

```rust
// Before:
pub is_contig_end: bool,

// After:
pub contig_end_genomes: HashSet<GenomeId>,
```

Change the `from_cluster` initialization:

```rust
// Before:
is_contig_end: false,

// After:
contig_end_genomes: HashSet::new(),
```

Same for `from_cluster_with_genes`.

- [ ] **Step 2: Update GraphBuilder to populate per-genome contig-end data**

In `src/graph/builder.rs`, change the contig-end marking logic. Currently at ~line 133-138:

```rust
// Before:
// Set is_contig_end if any gene in this cluster is at a contig boundary
node.is_contig_end = true;
break;  // stops at first contig-end gene

// After:
// Track which genomes have this cluster at a contig boundary
node.contig_end_genomes.insert(genome_id.clone());
// Don't break — collect ALL genomes with contig-end genes
```

- [ ] **Step 3: Update ContigEndPruner in contig_end.rs**

Change `src/correction/contig_end.rs:49`:

```rust
// Before:
node.is_contig_end && graph.is_degree_one(entry.key()) && node.support < self.min_support

// After:
!node.contig_end_genomes.is_empty() && graph.is_degree_one(entry.key()) && node.support < self.min_support
```

Update all test fixtures that set `is_contig_end = true`:

```rust
// Before:
end_node.is_contig_end = true;

// After:
end_node.contig_end_genomes.insert(GenomeId::new("genome1"));
```

Repeat for lines 114, 150, 187.

- [ ] **Step 4: Update MisassemblyEdgeCleaner in misassembly.rs**

Change `src/correction/misassembly.rs:67`:

```rust
// Before:
if node.is_contig_end {

// After:
if !node.contig_end_genomes.is_empty() {
```

Update test at line 137:

```rust
// Before:
end_node.is_contig_end = true;

// After:
end_node.contig_end_genomes.insert(GenomeId::new("genome1"));
```

- [ ] **Step 5: Update GraphBuilder tests in builder.rs**

Change `src/graph/builder.rs:336-338`:

```rust
// Before:
assert!(node_c1.is_contig_end, ...)
assert!(node_c2.is_contig_end, ...)
assert!(node_c3.is_contig_end, ...)

// After:
assert!(node_c1.contig_end_genomes.contains(&GenomeId::new("genome1")), ...)
assert!(node_c2.contig_end_genomes.contains(&GenomeId::new("genome1")), ...)
assert!(node_c3.contig_end_genomes.contains(&GenomeId::new("genome1")), ...)
```

- [ ] **Step 6: Update merge_nodes in concurrent.rs (already has contig_end_genomes from Task 2)**

The `merge_nodes` method already extends `contig_end_genomes` from Task 2. Verify this line is present:

```rust
target_node.contig_end_genomes.extend(source_node.contig_end_genomes);
```

- [ ] **Step 7: Update GML writer to serialize contig_end_genomes**

In `src/output/graph.rs`, change the GML serialization of the old `is_contig_end` attribute:

```rust
// Before (if present):
// is_contig_end true/false

// After:
// contig_end_genomes ["genome1","genome2",...]
```

Serialize as a comma-separated string of genome IDs.

- [ ] **Step 8: Run all tests**

Run: `cargo test --features full 2>&1 | head -80`
Expected: All tests pass with `contig_end_genomes: HashSet<GenomeId>`.

- [ ] **Step 9: Commit is_contig_end changes**

```bash
git add src/graph/types.rs src/graph/builder.rs src/graph/concurrent.rs src/correction/contig_end.rs src/correction/misassembly.rs src/output/graph.rs
git commit -m "refactor: change Node.is_contig_end to contig_end_genomes HashSet<GenomeId>

Per-gene contig-end tracking: nodes now track which specific genomes
have this cluster at a contig boundary, enabling more precise correction."
```

---

### Task 4: Add GML backward compatibility for old format

**Files:**
- Modify: `src/output/graph.rs` (GML reader)

- [ ] **Step 1: Add backward compatibility in GML reader**

When reading GML files, check for the old `is_contig_end` boolean attribute and the old `centroid_sequence` string attribute. If found, convert:

```rust
// Old format: is_contig_end true/false
// Convert: if is_contig_end is true, set contig_end_genomes to all member genomes
// Old format: centroid_sequence "ATCG..."
// Convert: wrap into centroid_sequences: vec![seq]
```

- [ ] **Step 2: Test GML round-trip with old format**

Write a test that reads a GML file with `is_contig_end true` and `centroid_sequence "ATCG"`, verifying that the resulting Node has `contig_end_genomes` populated and `centroid_sequences` as a `vec![]`.

- [ ] **Step 3: Commit backward compatibility**

```bash
git add src/output/graph.rs
git commit -m "feat: add GML backward compatibility for is_contig_end and centroid_sequence

Reads old GML format (boolean is_contig_end, single centroid_sequence) and
converts to new format (contig_end_genomes HashSet, centroid_sequences Vec)."
```