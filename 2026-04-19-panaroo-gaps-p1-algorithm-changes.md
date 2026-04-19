# Panaroo Feature Parity — Phase 1: Core Algorithm Changes

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement iterative multi-threshold collapsing, length-based filtering, consensus removal in refinding, and shortest-path paralog resolution.

**Architecture:** These are core algorithm changes that improve PanMiner's accuracy to match Panaroo. Iterative collapsing replaces single-pass at 0.70 with multi-threshold [0.99, 0.95, 0.9, 0.8, 0.7]. Length filtering adds a `len_dif_percent` gate in both clusterers. Consensus removal adds a `remove_by_consensus` mode to MissingGeneRecoverer. Shortest path adds `petgraph`-style BFS as primary paralog resolution.

**Tech Stack:** Rust, existing crate dependencies

---

### Task 1: Add iterative multi-threshold collapsing to FragmentMerger and pipeline

**Files:**
- Modify: `src/config.rs` (add `collapse_thresholds: Vec<f32>`)
- Modify: `src/correction/fragment.rs` (iterative collapsing loop)
- Modify: `src/pipeline.rs` (wire iterative collapsing)

- [ ] **Step 1: Add `collapse_thresholds` to PanminerConfig**

In `src/config.rs`, add the field and CLI argument:

```rust
/// Iterative collapsing thresholds (high to low).
/// Defaults to [0.99, 0.95, 0.9, 0.8, 0.7] matching Panaroo's progressive approach.
pub collapse_thresholds: Vec<f32>,
```

Add to CLI in `src/main.rs`:

```rust
/// Collapsing thresholds (comma-separated, high to low)
#[arg(long, default_value = "0.99,0.95,0.9,0.8,0.7")]
collapse_thresholds: String,
```

Parse in the CLI handler:

```rust
let collapse_thresholds: Vec<f32> = cli.collapse_thresholds
    .split(',')
    .filter_map(|s| s.parse().ok())
    .collect();
```

- [ ] **Step 2: Modify FragmentMerger to support configurable threshold**

In `src/correction/fragment.rs`, `FragmentMerger` already has `collapse_threshold: f32`. Replace it with `collapse_thresholds: Vec<f32>`:

```rust
pub struct FragmentMerger {
    coverage_threshold: f32,
    identity_threshold: f32,
    collapse_thresholds: Vec<f32>,  // NEW: replaces single collapse_threshold
    bfs_depth: usize,
}
```

Update constructors:

```rust
impl FragmentMerger {
    pub fn new() -> Self {
        Self {
            coverage_threshold: 0.95,
            identity_threshold: 0.99,
            collapse_thresholds: vec![0.99, 0.95, 0.9, 0.8, 0.7],
            bfs_depth: 3,
        }
    }

    pub fn with_collapse_thresholds(mut self, thresholds: Vec<f32>) -> Self {
        self.collapse_thresholds = thresholds;
        self
    }
}
```

Keep `with_collapse_threshold` as a compatibility shim that sets a single threshold:

```rust
pub fn with_collapse_threshold(mut self, threshold: f32) -> Self {
    self.collapse_thresholds = vec![threshold];
    self
}
```

- [ ] **Step 3: Modify `collapse_gene_families_with_cache` to accept a threshold parameter**

Change the method signature to accept an explicit threshold:

```rust
pub fn collapse_gene_families_with_threshold(
    &self,
    graph: &ConcurrentGraph,
    sequences: &HashMap<String, Vec<u8>>,
    threshold: f32,
    cache: Option<&mut DistanceCache>,
) -> Result<usize>
```

The method body uses `threshold` instead of `self.collapse_threshold`. Keep `collapse_gene_families_with_cache` as a wrapper that calls this with `self.collapse_thresholds[0]`.

- [ ] **Step 4: Update pipeline.rs to use iterative collapsing**

In `src/pipeline.rs`, `run_corrections`, replace the two single-collapse calls with an iterative loop:

```rust
// Phase 4.3: Mistranslation correction (stays the same)
let merger = FragmentMerger::new()
    .with_collapse_thresholds(self.config.collapse_thresholds.clone());
let mistrans_stats = merger.correct_mistranslations(&graph, &sequences)?;

// Phase 4.4-4.6: Iterative gene family collapsing
let mut distance_cache = DistanceCache::new();
let mut total_merged = 0usize;
for threshold in &merger.collapse_thresholds {
    let merged = merger.collapse_gene_families_with_threshold(
        &graph, &sequences, *threshold, Some(&mut distance_cache)
    )?;
    total_merged += merged;
    if merged == 0 {
        break; // No more merges possible at this threshold
    }
}

// Phase 4.5: Missing gene recovery
let recoverer = MissingGeneRecoverer::new()
    .with_min_identity(0.70)
    .with_search_window(5000);
let recovery_stats = recoverer.recover(&graph, &contig_sequences, &cluster_sequences)?;

// Phase 4.6: Re-collapse families after recovery
let merger2 = FragmentMerger::new()
    .with_collapse_thresholds(self.config.collapse_thresholds.clone());
for threshold in &merger2.collapse_thresholds {
    let merged = merger2.collapse_gene_families_with_threshold(
        &graph, &sequences_after_recovery, *threshold, Some(&mut distance_cache)
    )?;
    if merged == 0 {
        break;
    }
}
```

- [ ] **Step 5: Write tests for iterative collapsing**

Add a test in `src/correction/fragment.rs` that creates a graph with clusters at different identity levels and verifies that multi-threshold collapsing merges more clusters than single-threshold:

```rust
#[test]
fn test_iterative_collapsing_merges_more() {
    // Create clusters at 85%, 75%, and 60% identity
    // Single threshold at 0.70 should merge only the 75% pair
    // Iterative [0.99, 0.95, 0.9, 0.8, 0.7] should merge all three
    // ...
}
```

- [ ] **Step 6: Run all tests**

Run: `cargo test --features full 2>&1 | head -80`
Expected: All tests pass including new iterative collapsing test.

- [ ] **Step 7: Commit iterative collapsing**

```bash
git add src/config.rs src/correction/fragment.rs src/pipeline.rs src/main.rs
git commit -m "feat: iterative multi-threshold gene family collapsing

Replace single-pass collapse at 0.70 with progressive collapsing at
[0.99, 0.95, 0.9, 0.8, 0.7] (configurable). Matches Panaroo's collapse_families
behavior for improved distant homolog detection. DistanceCache persists across
thresholds for efficiency."
```

---

### Task 2: Add length-based filtering during clustering

**Files:**
- Modify: `src/config.rs` (add `len_dif_percent: f32`)
- Modify: `src/clustering/cpu.rs` (add length filter to `sequence_identity`)
- Modify: `src/clustering/mmseqs.rs` (add `--cov-mode` and `-c` flags)
- Modify: `src/main.rs` (add CLI argument)

- [ ] **Step 1: Add `len_dif_percent` to PanminerConfig**

In `src/config.rs`:

```rust
/// Length difference cutoff for clustering (0.0-1.0, default 0.98).
/// Gene pairs with length difference > (1 - len_dif_percent) are excluded from clustering.
/// Matches CD-HIT's -s parameter.
pub len_dif_percent: f32,
```

Default value: `0.98`.

Add CLI argument in `src/main.rs`:

```rust
/// Length difference cutoff for clustering (0.0-1.0)
#[arg(long, default_value = "0.98")]
len_dif_percent: f32,
```

- [ ] **Step 2: Add length filtering to CpuClusterer**

In `src/clustering/cpu.rs`, modify the `greedy_cluster` method to add a length filter before identity comparison:

```rust
fn greedy_cluster(&self, genes: &[Gene], identity_threshold: f32, len_dif_percent: f32) -> Vec<GeneCluster> {
    // ... for each gene:
    // Find best matching centroid
    for (centroid_idx, centroid) in centroids.iter().enumerate() {
        // Length filter: skip if length difference exceeds threshold
        let max_len = centroid.sequence.len().max(gene.sequence.len()) as f32;
        let len_diff = (centroid.sequence.len().abs_diff(gene.sequence.len())) as f32 / max_len;
        if len_diff > (1.0 - len_dif_percent) {
            continue; // Skip this pair
        }

        // Identity check (existing logic)
        let identity = Self::sequence_identity(&gene.sequence, &centroid.sequence);
        if identity >= identity_threshold {
            // Add to cluster
            // ...
        }
    }
}
```

Update the `Clusterer` trait to pass `len_dif_percent`:

```rust
pub trait Clusterer {
    fn cluster(&self, genes: &[Gene], identity_threshold: f32, len_dif_percent: f32) -> Result<Vec<GeneCluster>>;
    fn name(&self) -> &str;
    fn is_available(&self) -> bool;
}
```

- [ ] **Step 3: Add length filtering flags to MMseqs2**

In `src/clustering/mmseqs.rs`, update `easy_cluster` to pass coverage flags:

```rust
// In the Command construction:
.arg("--cov-mode").arg("1")       // Coverage of shorter sequence
.arg("-c").arg(len_dif_percent.to_string())  // Min coverage
```

Update `MMseqsRunner::cluster` signature to accept `len_dif_percent`.

- [ ] **Step 4: Update pipeline to pass `len_dif_percent`**

In `src/pipeline.rs`, update `cluster_genes` to pass `self.config.len_dif_percent` to both clusterers.

- [ ] **Step 5: Write test for length filtering**

```rust
#[test]
fn test_length_filter_rejects_different_lengths() {
    // Create two genes: one 300bp, one 150bp
    // At len_dif_percent=0.98, max allowed difference is 2%
    // 150/300 = 50% difference, should be rejected
    // At len_dif_percent=0.50, should be accepted
}
```

- [ ] **Step 6: Run all tests**

Run: `cargo test --features full 2>&1 | head -80`
Expected: All tests pass.

- [ ] **Step 7: Commit length filtering**

```bash
git add src/config.rs src/clustering/cpu.rs src/clustering/mmseqs.rs src/clustering/traits.rs src/pipeline.rs src/main.rs
git commit -m "feat: add length-based filtering during clustering

Gene pairs with length difference > (1 - len_dif_percent) are now excluded
from clustering, matching CD-HIT's -s parameter behavior. Default 0.98 means
sequences differing by >2% in length are separated. MMseqs2 uses --cov-mode 1
and -c flags for the same effect."
```

---

### Task 3: Add consensus removal to MissingGeneRecoverer

**Files:**
- Modify: `src/correction/missing.rs` (add `remove_by_consensus` field and logic)
- Modify: `src/config.rs` (add mode-dependent default)
- Modify: `src/pipeline.rs` (wire consensus removal)

- [ ] **Step 1: Add `remove_by_consensus` field to MissingGeneRecoverer**

In `src/correction/missing.rs`:

```rust
pub struct MissingGeneRecoverer {
    min_identity: f32,
    search_window: usize,
    prop_match: f32,
    remove_by_consensus: bool,  // NEW: delete nodes where refound hits exceed original size
}
```

Update constructors:

```rust
impl MissingGeneRecoverer {
    pub fn new() -> Self {
        Self {
            min_identity: 0.70,
            search_window: 5000,
            prop_match: 0.20,
            remove_by_consensus: false,
        }
    }

    pub fn with_remove_by_consensus(mut self, remove: bool) -> Self {
        self.remove_by_consensus = remove;
        self
    }
}
```

- [ ] **Step 2: Implement consensus removal logic**

In the `recover` method, after adding refound genes, add:

```rust
if self.remove_by_consensus {
    // For each node that received refound hits:
    // If total_refound_hits > original_node.support, mark node for removal
    let nodes_to_remove: Vec<ClusterId> = recovery_counts
        .iter()
        .filter(|(id, (original_size, refound_count))| refound_count > original_size)
        .map(|(id, _)| id.clone())
        .collect();

    for node_id in nodes_to_remove {
        graph.remove_node(&node_id);
    }
}
```

- [ ] **Step 3: Wire consensus removal into pipeline based on correction mode**

In `src/pipeline.rs`, `run_corrections`:

```rust
let remove_by_consensus = matches!(self.config.correction_mode, CorrectionMode::Strict);
let recoverer = MissingGeneRecoverer::new()
    .with_min_identity(0.70)
    .with_search_window(5000)
    .with_remove_by_consensus(remove_by_consensus);
```

- [ ] **Step 4: Write test for consensus removal**

```rust
#[test]
fn test_consensus_removes_spurious_nodes() {
    // Create a graph where a node with support 2 gets 5 refound hits
    // With remove_by_consensus=true, the node should be removed
    // With remove_by_consensus=false, the node should remain
}
```

- [ ] **Step 5: Run all tests**

Run: `cargo test --features full 2>&1 | head -80`
Expected: All tests pass.

- [ ] **Step 6: Commit consensus removal**

```bash
git add src/correction/missing.rs src/config.rs src/pipeline.rs
git commit -m "feat: add consensus removal in missing gene recovery

When remove_by_consensus is enabled (strict mode), nodes where refound hits
exceed the original node size are deleted. Matches Panaroo's --remove_by_consensus
behavior for eliminating spurious annotation artifacts."
```

---

### Task 4: Add shortest-path paralog resolution

**Files:**
- Modify: `src/correction/paralog.rs` (add BFS shortest path as primary method)

- [ ] **Step 1: Implement shortest path distance function**

In `src/correction/paralog.rs`, add a function that computes shortest path distance between two nodes using BFS on the ConcurrentGraph:

```rust
/// Compute shortest path distance between two nodes using BFS.
/// Returns None if no path exists within max_depth.
fn shortest_path_distance(
    graph: &ConcurrentGraph,
    from: &ClusterId,
    to: &ClusterId,
    max_depth: usize,
) -> Option<usize> {
    if from == to {
        return Some(0);
    }
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();
    visited.insert(from.clone());
    queue.push_back((from.clone(), 0usize));

    while let Some((current, depth)) = queue.pop_front() {
        if depth >= max_depth {
            continue;
        }
        for neighbor in graph.neighbors(&current) {
            if neighbor == *to {
                return Some(depth + 1);
            }
            if visited.insert(neighbor.clone()) {
                queue.push_back((neighbor, depth + 1));
            }
        }
    }
    None
}
```

- [ ] **Step 2: Integrate shortest path into paralog resolution**

In the `resolve` method of `ParalogResolver`, use shortest path as the primary method:

```rust
// For each paralog copy, try shortest path first
for (para_id, reference_id) in &assignments {
    if let Some(distance) = shortest_path_distance(graph, para_id, reference_id, self.max_context) {
        // Use shortest path distance for assignment
        assignment_scores.insert(para_id.clone(), (reference_id.clone(), 1.0 / (1.0 + distance as f64)));
    } else {
        // Fall back to context vector similarity
        let context_sim = compute_context_similarity(graph, para_id, reference_id, self.max_context);
        assignment_scores.insert(para_id.clone(), (reference_id.clone(), context_sim));
    }
}
```

- [ ] **Step 3: Write test for shortest path resolution**

```rust
#[test]
fn test_shortest_path_resolves_paralogs() {
    // Create a graph with two paralog nodes connected by a short path
    // Verify that shortest_path_distance returns the correct distance
    // Verify that context vector similarity is used as fallback when no path exists
}
```

- [ ] **Step 4: Run all tests**

Run: `cargo test --features full 2>&1 | head -80`
Expected: All tests pass.

- [ ] **Step 5: Commit shortest path paralog resolution**

```bash
git add src/correction/paralog.rs
git commit -m "feat: add shortest path as primary paralog resolution method

BFS shortest path distance is now the primary method for resolving paralog
copies, matching Panaroo's nx.shortest_path_length approach. Context vector
similarity (BFS depth 5) is used as fallback when no path exists within
max_context depth."
```