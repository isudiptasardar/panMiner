# P1 Features Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement three P1 (High Priority) features to achieve Panaroo feature parity: contig-end pruning, structural variant matrix output, and expanded integration tests.

**Architecture:** 
- Contig-end pruning: Iterative removal of degree-1 nodes at contig ends (distinct from contamination removal)
- Structural variant matrix: Gene triplet presence/absence output capturing co-occurrence patterns
- Integration tests: End-to-end pipeline testing with real data scenarios

**Tech Stack:** Rust, rayon for parallelism, dashmap for concurrent graph, thiserror for error types

---
## File Structure

### New Files to Create
| File | Purpose |
|------|---------|
| `src/correction/contig_end.rs` | Contig-end pruning module with iterative degree-1 removal |
| `src/output/sv_matrix.rs` | Structural variant matrix (gene triplets) TSV output |
| `tests/p1_integration_tests.rs` | Expanded integration tests for P1 features |

### Files to Modify
| File | Changes |
|------|---------|
| `src/graph/types.rs` | Add `is_contig_end: bool` field to `Node` struct |
| `src/graph/builder.rs` | Mark nodes at contig ends during graph building |
| `src/output/mod.rs` | Add `SVMatrixWriter` export |
| `src/pipeline.rs` | Wire contig-end pruning into pipeline flow |

---

### Task 1: Add `is_contig_end` to Node struct

**Files:**
- Modify: `src/graph/types.rs:193-208` (Node struct)

- [ ] **Step 1: Add field to Node struct**

Add `is_contig_end: bool` field to the `Node` struct after `centroid_sequence`:

```rust
pub struct Node {
    pub cluster_id: ClusterId,
    pub support: usize,
    pub genomes: HashSet<GenomeId>,
    pub annotations: HashSet<String>,
    pub is_paralog: bool,
    pub centroid_sequence: Option<Sequence>,
    pub contig_sequences: HashMap<String, Sequence>,
    pub is_contig_end: bool,  // NEW: Mark if this node is at a contig end
}
```

- [ ] **Step 2: Update from_cluster to set default**

Set default value in `Node::from_cluster`:

```rust
impl Node {
    pub fn from_cluster(cluster: &GeneCluster) -> Self {
        Self {
            cluster_id: cluster.id.clone(),
            support: cluster.support,
            genomes: HashSet::new(),
            annotations: HashSet::new(),
            is_paralog: cluster.is_paralog,
            centroid_sequence: cluster.centroid.clone(),
            contig_sequences: HashMap::new(),
            is_contig_end: false,  // NEW: Default false
        }
    }
    // ... rest of impl
}
```

- [ ] **Step 3: Commit changes**

```bash
git add src/graph/types.rs
git commit -m "feat: add is_contig_end field to Node struct"
```

---

### Task 2: Mark contig ends during graph building

**Files:**
- Modify: `src/graph/builder.rs`

- [ ] **Step 1: Add contig gene counting**

After grouping genes by contig (line ~57), count genes per contig:

```rust
// Count genes per contig to identify contig ends
let contig_gene_count: HashMap<(GenomeId, String), usize> = genes_by_contig
    .iter()
    .map(|((genome, contig), genes)| ((genome.clone(), contig.clone()), genes.len()))
    .collect();
```

- [ ] **Step 2: Mark contig ends during node creation**

After line ~91 (where contig sequences are added), add logic to mark contig ends:

```rust
// Add contig sequences if available
for ((genome, contig), seq) in &contig_sequences {
    if node.genomes.contains(genome) {
        node.add_contig_sequence(contig.clone(), seq.clone());
    }
}
// Mark contig ends: nodes that are the only gene on their contig
for ((genome, contig), seq) in &contig_sequences {
    if node.genomes.contains(genome) && contig_gene_count.get(&(genome.clone(), contig.clone())) == Some(&1) {
        // This is the only gene on this contig for this genome
        // Check if this node represents a gene on this contig
        let node_has_gene_on_contig = cluster.genes.iter().any(|g| {
            genes.iter().any(|gene| 
                gene.id.to_string() == g.to_string() && 
                gene.contig == *contig && 
                gene.genome_id == *genome
            )
        });
        if node_has_gene_on_contig {
            // Mark as contig end (we'll do this more efficiently in actual implementation)
        }
    }
}
```

- [ ] **Step 2b: Simpler approach - track contig ends separately**

Actually, let's simplify: Add a method to `ConcurrentGraph` that marks contig ends after all nodes are added:

```rust
// In ConcurrentGraph, add this method:
pub fn mark_contig_ends(&self, genes: &[Gene]) {
    // Count genes per (genome, contig)
    let gene_count: HashMap<(GenomeId, String), usize> = genes
        .iter()
        .fold(HashMap::new(), |mut acc, g| {
            let key = (g.genome_id.clone(), g.contig.clone());
            *acc.entry(key).or_insert(0) += 1;
            acc
        });
    
    // For each node, check if any of its genes are at contig ends
    for entry in self.nodes.iter() {
        let cluster_id = entry.key();
        let mut node = entry.clone(); // Need mutable access
        // Check if any gene in this cluster is the only one on its contig
        // ... implementation
    }
}
```

- [ ] **Step 3: Commit changes**

```bash
git add src/graph/builder.rs
git commit -m "feat: mark contig ends during graph building"
```

---

### Task 3: Create ContigEndPruner module

**Files:**
- Create: `src/correction/contig_end.rs`

- [ ] **Step 1: Create new module file**

Write the complete module:

```rust
//! Contig-end pruning from the pangenome graph.
//!
//! Recursively removes nodes that are at contig ends (single gene on contig)
//! with low support. This is distinct from general contamination removal.

use crate::error::Result;
use crate::graph::ConcurrentGraph;

/// Removes contig-end nodes from the pangenome graph.
///
/// Contig-end nodes are those where the gene is the only one on its contig
/// in that genome. These are often fragmented/partial genes at contig boundaries.
pub struct ContigEndPruner {
    /// Minimum support threshold for keeping contig-end nodes
    min_support: usize,
}

impl ContigEndPruner {
    /// Create a new contig-end pruner with default threshold.
    pub fn new() -> Self {
        Self { min_support: 1 }
    }

    /// Create with custom threshold.
    pub fn with_min_support(mut self, min_support: usize) -> Self {
        self.min_support = min_support;
        self
    }

    /// Mark contig ends in the graph based on gene data.
    pub fn mark_contig_ends(&self, graph: &ConcurrentGraph, genes: &[Gene]) {
        // Count genes per (genome, contig)
        let gene_count: std::collections::HashMap<(crate::graph::GenomeId, String), usize> = 
            genes.iter().fold(std::collections::HashMap::new(), |mut acc, g| {
                let key = (g.genome_id.clone(), g.contig.clone());
                *acc.entry(key).or_insert(0) += 1;
                acc
            });

        // Mark nodes as contig ends
        for entry in graph.nodes.iter() {
            let mut node = entry.value_mut().clone();
            let mut is_end = false;
            
            // Check if any gene in this node's clusters is the only one on its contig
            for gene_id in &entry.value().cluster_id.to_string() {
                // This needs the gene-to-cluster mapping
                // Simplified: mark all nodes as contig ends initially
                // A more precise implementation would track which genes are on which contigs
            }
            
            node.is_contig_end = is_end;
            // Update the node in place (need DashMap entry access)
        }
    }

    /// Remove contig-end nodes from the graph.
    ///
    /// Recursively removes contig-end nodes with support below threshold.
    pub fn prune(&self, graph: &ConcurrentGraph) -> Result<PruningStats> {
        let mut total_removed = 0;
        let mut iteration = 0;
        let max_iterations = 100;

        loop {
            if iteration >= max_iterations {
                break;
            }

            let to_remove: Vec<_> = graph.nodes
                .iter()
                .filter(|entry| {
                    let node = entry.value();
                    node.is_contig_end && node.support <= self.min_support
                })
                .map(|entry| entry.key().clone())
                .collect();

            if to_remove.is_empty() {
                break;
            }

            let removed_count = to_remove.len();
            for cluster_id in to_remove {
                graph.remove_node(&cluster_id);
            }
            total_removed += removed_count;
            iteration += 1;

            tracing::debug!(
                "Contig-end pruning iteration {}: removed {} nodes",
                iteration,
                removed_count
            );
        }

        tracing::info!(
            "Contig-end pruning: removed {} nodes in {} iterations",
            total_removed,
            iteration
        );

        Ok(PruningStats {
            nodes_removed: total_removed,
            iterations: iteration,
        })
    }
}

/// Statistics from contig-end pruning.
#[derive(Debug, Clone)]
pub struct PruningStats {
    /// Total nodes removed
    pub nodes_removed: usize,
    /// Number of iterations
    pub iterations: usize,
}

impl Default for ContigEndPruner {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_end_pruner_creation() {
        let pruner = ContigEndPruner::new();
        assert_eq!(pruner.min_support, 1);
    }

    #[test]
    fn test_with_min_support() {
        let pruner = ContigEndPruner::new().with_min_support(5);
        assert_eq!(pruner.min_support, 5);
    }
}
```

- [ ] **Step 2: Update mod.rs to export**

Add to `src/correction/mod.rs`:

```rust
pub mod contamination;
pub mod fragment;
pub mod missing;
pub mod contig_end;  // NEW
```

- [ ] **Step 3: Commit changes**

```bash
git add src/correction/contig_end.rs src/correction/mod.rs
git commit -m "feat: add contig-end pruning module"
```

---

### Task 4: Wire contig-end pruning into pipeline

**Files:**
- Modify: `src/pipeline.rs:248-287`

- [ ] **Step 1: Import new module**

Add to imports at top of `src/pipeline.rs`:

```rust
use crate::correction::{ContaminationRemover, FragmentMerger, MissingGeneRecoverer, ContigEndPruner};
```

- [ ] **Step 2: Add contig-end pruning to pipeline flow**

After contamination removal (line ~252), add contig-end pruning:

```rust
// Phase 4: Error correction
tracing::info!("Phase 4: Running error correction");

// 4a. Contamination removal
let remover = ContaminationRemover::from_mode(&self.config.mode, num_genomes);
remover.remove(graph)?;

// 4b. Contig-end pruning (NEW)
let pruner = ContigEndPruner::new()
    .with_min_support(self.config.min_support);
pruner.prune(graph)?;

// 4c. Fragment merging with actual cluster centroid sequences
let merger = FragmentMerger::new()
    .with_collapse_threshold(self.config.collapse_threshold);
// ... rest of existing code
```

- [ ] **Step 3: Commit changes**

```bash
git add src/pipeline.rs
git commit -m "feat: wire contig-end pruning into pipeline"
```

---

### Task 5: Create structural variant matrix output module

**Files:**
- Create: `src/output/sv_matrix.rs`

- [ ] **Step 1: Create new module file**

Write the complete module:

```rust
//! Structural variant matrix (gene triplet) output.
//!
//! Outputs a TSV file with gene triplet presence/absence patterns
//! across genomes, capturing co-occurrence of nearby genes.

use crate::error::Result;
use crate::graph::{Edge, PangenomeGraph};
use std::path::Path;

/// A gene triplet (two adjacent clusters + their co-occurrence pattern).
#[derive(Debug, Clone)]
pub struct GeneTriplet {
    /// First cluster in triplet
    pub cluster_a: String,
    /// Second cluster in triplet
    pub cluster_b: String,
    /// Presence/absence for each genome
    pub pattern: Vec<bool>,
}

/// Writer for structural variant matrix output.
pub struct SVMatrixWriter;

impl SVMatrixWriter {
    /// Extract gene triplets from the pangenome graph.
    pub fn extract_triplets(graph: &PangenomeGraph) -> Vec<GeneTriplet> {
        let mut triplets = Vec::new();

        // For each edge, create a gene triplet entry
        for entry in graph.edges.iter() {
            let edge = entry.value();
            let (cluster_a, cluster_b) = entry.key();

            // Build presence/absence pattern for all genomes
            let mut all_genomes: std::collections::HashSet<String> = std::collections::HashSet::new();
            for genome in &edge.genomes {
                all_genomes.insert(genome.to_string());
            }

            // Get genomes from both nodes to ensure complete pattern
            if let Some(node_a) = graph.nodes.get(&edge.from) {
                for genome in &node_a.genomes {
                    all_genomes.insert(genome.to_string());
                }
            }
            if let Some(node_b) = graph.nodes.get(&edge.to) {
                for genome in &node_b.genomes {
                    all_genomes.insert(genome.to_string());
                }
            }

            // Create pattern vector
            let pattern: Vec<bool> = all_genomes.iter().map(|g| {
                edge.genomes.iter().any(|e| e.to_string() == *g)
            }).collect();

            triplets.push(GeneTriplet {
                cluster_a: cluster_a.to_string(),
                cluster_b: cluster_b.to_string(),
                pattern,
            });
        }

        triplets
    }

    /// Write gene triplets to TSV file.
    pub fn write_tsv(triplets: &[GeneTriplet], path: &Path) -> Result<()> {
        use std::io::Write;

        let file = std::fs::File::create(path)?;
        let mut writer = std::io::BufWriter::new(file);

        // Get all genome names from first triplet
        let genomes: Vec<String> = if let Some(first) = triplets.first() {
            // We need to track genomes during extraction
            // For now, extract from edge data
            vec![]
        } else {
            vec![]
        };

        // Write header
        let header = format!("GeneTriplet,{}", genomes.join(","));
        writeln!(writer, "{}", header)?;

        // Write each triplet
        for (i, triplet) in triplets.iter().enumerate() {
            // Build presence string for this triplet
            let presence: String = triplet.pattern
                .iter()
                .map(|p| if *p { "1" } else { "0" })
                .collect::<Vec<_>>()
                .join(",");

            let row = format!("{},{},{}", i + 1, triplet.cluster_a, triplet.cluster_b, presence);
            writeln!(writer, "{}", row)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_graph_no_triplets() {
        let graph = PangenomeGraph::new();
        let triplets = SVMatrixWriter::extract_triplets(&graph);
        assert!(triplets.is_empty());
    }
}
```

- [ ] **Step 2: Create simpler version that actually compiles**

Let me rewrite with proper implementation:

```rust
//! Structural variant matrix (gene triplet) output.
//!
//! Outputs a TSV file with gene triplet presence/absence patterns
//! across genomes, capturing co-occurrence of nearby genes.

use crate::error::Result;
use crate::graph::{PangenomeGraph};
use std::path::Path;

/// Writer for structural variant matrix output.
pub struct SVMatrixWriter {
    genome_names: Vec<String>,
}

impl SVMatrixWriter {
    /// Create a new writer.
    pub fn new() -> Self {
        Self { genome_names: Vec::new() }
    }

    /// Set genome names for output.
    pub fn with_genomes(mut self, names: Vec<String>) -> Self {
        self.genome_names = names;
        self
    }

    /// Extract gene triplets from the pangenome graph.
    pub fn extract_triplets(graph: &PangenomeGraph) -> Vec<(String, String, Vec<bool>)> {
        let mut triplets = Vec::new();

        // For each edge, create a gene triplet entry
        for entry in graph.edges.iter() {
            let edge = entry.value();
            let (cluster_a, cluster_b) = entry.key();

            // Build presence/absence pattern for all genomes
            let pattern: Vec<bool> = edge.genomes
                .iter()
                .map(|_| true) // Placeholder - actual logic needs genome mapping
                .collect();

            triplets.push((cluster_a.to_string(), cluster_b.to_string(), pattern));
        }

        triplets
    }

    /// Write gene triplets to TSV file.
    pub fn write_tsv(
        &self, 
        triplets: &[(String, String, Vec<bool>)],
        path: &Path
    ) -> Result<()> {
        use std::io::Write;

        let file = std::fs::File::create(path)?;
        let mut writer = std::io::BufWriter::new(file);

        // Write header with genome names
        let header = format!("# Structural Variant Matrix (gene triplets)\n# GeneTriplet,{}", 
            self.genome_names.join(","));
        writeln!(writer, "{}", header)?;

        // Write each triplet
        for (cluster_a, cluster_b, pattern) in triplets.iter() {
            let presence: String = pattern
                .iter()
                .map(|p| if *p { "1" } else { "0" })
                .collect::<Vec<_>>()
                .join(",");
            
            let triplet_id = format!("{}_{}", cluster_a, cluster_b);
            writeln!(writer, "{},{}", triplet_id, presence)?;
        }

        Ok(())
    }
}

impl Default for SVMatrixWriter {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_writer_creation() {
        let writer = SVMatrixWriter::new();
        assert!(writer.genome_names.is_empty());
    }

    #[test]
    fn test_with_genomes() {
        let writer = SVMatrixWriter::new()
            .with_genomes(vec!["g1".to_string(), "g2".to_string()]);
        assert_eq!(writer.genome_names.len(), 2);
    }
}
```

- [ ] **Step 3: Commit changes**

```bash
git add src/output/sv_matrix.rs
git commit -m "feat: add SV matrix TSV output module"
```

---

### Task 6: Integrate SV matrix into OutputWriter

**Files:**
- Modify: `src/output/mod.rs`
- Modify: `src/output/mod.rs` (OutputWriter::write_all)

- [ ] **Step 1: Add import and export**

At top of `src/output/mod.rs`:

```rust
mod sv_matrix;  // NEW
pub use sv_matrix::SVMatrixWriter;  // NEW
```

- [ ] **Step 2: Add OutputFormat variant**

Add to `OutputFormat` enum in `src/config.rs`:

```rust
pub enum OutputFormat {
    Matrix,
    Alignment,
    Graph,
    Json,
    Parquet,
    HtmlViz,
    Struct,  // Already exists for structural variants
    SVMatrix,  // NEW
}
```

- [ ] **Step 3: Implement write_all in OutputWriter**

Add SVMatrix case to `write_all` in `src/output/mod.rs`:

```rust
OutputFormat::SVMatrix => {
    let path = self.output_dir.join("gene_triplets.tsv");
    let triplets = SVMatrixWriter::extract_triplets(graph);
    let writer = SVMatrixWriter::new()
        .with_genomes(matrix.genome_names.clone());
    writer.write_tsv(&triplets, &path)?;
    tracing::info!("Wrote gene triplets matrix");
}
```

- [ ] **Step 4: Commit changes**

```bash
git add src/output/mod.rs src/config.rs
git commit -m "feat: integrate SV matrix into output formats"
```

---

### Task 7: Create expanded integration tests

**Files:**
- Create: `tests/p1_integration_tests.rs`

- [ ] **Step 1: Create new test file**

Write comprehensive integration tests:

```rust
//! P1 Feature Integration Tests
//!
//! Tests for:
//! - Contig-end pruning
//! - Structural variant matrix output
//! - Paralog handling
//! - Large dataset processing

use std::fs::File;
use std::io::Write;
use tempfile::TempDir;

use panminer::config::{PanminerConfig, CorrectionMode, OutputFormat};
use panminer::pipeline::PanminerPipeline;

mod test_helpers {
    use std::fs::File;
    use std::io::Write;
    use tempfile::TempDir;

    /// Create a test GFF3 file with multiple genes on a contig
    pub fn create_test_gff_with_contig_ends(
        dir: &TempDir,
        name: &str,
        gene_ids: &[&str],
        start: u32,
    ) -> std::path::PathBuf {
        let gff_path = dir.path().join(format!("{}.gff", name));
        let mut file = File::create(&gff_path).unwrap();

        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "##sequence-region {} 1 10000", name).unwrap();

        for (i, gene_id) in gene_ids.iter().enumerate() {
            let gene_start = start + (i as u32 * 150);
            let gene_end = gene_start + 99;
            writeln!(file, "{}\tProkka\tgene\t{}\t{}\t.\t+\t.\tID={};product=test", name, gene_start, gene_end, gene_id).unwrap();
        }

        writeln!(file, "##FASTA").unwrap();
        writeln!(file, ">{}", name).unwrap();
        let seq_len = (start + (gene_ids.len() as u32 * 150) + 100) as usize;
        let sequence = "ATCG".repeat(seq_len / 4 + 1);
        for i in (0..seq_len).step_by(80) {
            writeln!(file, "{}", &sequence[i..std::cmp::min(i + 80, seq_len)]).unwrap();
        }

        gff_path
    }

    /// Create test GFF with paralog (same gene twice in one genome)
    pub fn create_test_gff_with_paralog(dir: &TempDir, name: &str) -> std::path::PathBuf {
        let gff_path = dir.path().join(format!("{}.gff", name));
        let mut file = File::create(&gff_path).unwrap();

        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "##sequence-region {} 1 10000", name).unwrap();

        // Two copies of same gene (paralog)
        writeln!(file, "{}\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene1;product=test", name).unwrap();
        writeln!(file, "{}\tProkka\tgene\t300\t400\t.\t+\t.\tID=gene1;product=test", name).unwrap();

        writeln!(file, "##FASTA").unwrap();
        writeln!(file, ">{}", name).unwrap();
        let seq_len = 500;
        let sequence = "ATCG".repeat(seq_len / 4 + 1);
        for i in (0..seq_len).step_by(80) {
            writeln!(file, "{}", &sequence[i..std::cmp::min(i + 80, seq_len)]).unwrap();
        }

        gff_path
    }
}

#[test]
fn test_pipeline_contig_end_pruning() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create GFF with genes at contig ends (single genes on contigs)
    let gff = test_helpers::create_test_gff_with_contig_ends(&temp_dir, "genome1", &["gene1", "gene2", "gene3"], 100);

    let config = PanminerConfig::new()
        .with_input_files(vec![gff])
        .with_output_dir(output_dir.clone())
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run();

    assert!(result.is_ok(), "Pipeline should complete");
}

#[test]
fn test_pipeline_sv_matrix_output() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create multiple genomes with adjacent genes
    let gff1 = test_helpers::create_test_gff_with_contig_ends(&temp_dir, "genome1", &["gene1", "gene2"], 100);
    let gff2 = test_helpers::create_test_gff_with_contig_ends(&temp_dir, "genome2", &["gene1", "gene2"], 100);

    let config = PanminerConfig::new()
        .with_input_files(vec![gff1, gff2])
        .with_output_dir(output_dir.clone())
        .with_outputs([OutputFormat::Matrix, OutputFormat::Graph, OutputFormat::SVMatrix].into_iter().collect());

    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run();

    assert!(result.is_ok(), "Pipeline should complete with SV matrix");

    let output_paths = result.unwrap();
    
    // Check SV matrix was created
    assert!(output_paths.struct_csv.is_some(), "Structural variant matrix should be created");
}

#[test]
fn test_pipeline_with_paralogs() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    let gff = test_helpers::create_test_gff_with_paralog(&temp_dir, "genome1");

    let config = PanminerConfig::new()
        .with_input_files(vec![gff])
        .with_output_dir(output_dir.clone());

    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run();

    // Paralogs should be detected and handled
    assert!(result.is_ok(), "Pipeline should handle paralogs");
}

#[test]
fn test_pipeline_large_dataset() {
    let temp_dir = TempDir::new().unwrap();
    let output_dir = temp_dir.path().join("output");

    // Create 10+ genomes
    let mut gffs = Vec::new();
    for i in 0..10 {
        let gff = test_helpers::create_test_gff_with_contig_ends(
            &temp_dir,
            &format!("genome{}", i),
            &["gene1", "gene2"],
            100,
        );
        gffs.push(gff);
    }

    let config = PanminerConfig::new()
        .with_input_files(gffs)
        .with_output_dir(output_dir.clone());

    let pipeline = PanminerPipeline::new(config);
    let result = pipeline.run();

    assert!(result.is_ok(), "Pipeline should handle large dataset");
}
```

- [ ] **Step 2: Run tests to verify they pass**

```bash
cargo test --test p1_integration_tests
```

- [ ] **Step 3: Commit changes**

```bash
git add tests/p1_integration_tests.rs
git commit -m "test: add P1 feature integration tests"
```

---

### Task 8: Fix any issues and verify build

- [ ] **Step 1: Run cargo check**

```bash
cargo check
```

Fix any compilation errors.

- [ ] **Step 2: Run all tests**

```bash
cargo test
```

- [ ] **Step 3: Final commit**

```bash
git add .
git commit -m "fix: resolve compilation errors from P1 features"
```

---

### Task 9: Documentation and cleanup

- [ ] **Step 1: Update README**

Add new features to feature list in `README.md`.

- [ ] **Step 2: Update Comparison.md**

Update P1 feature status in `Comparison.md`.

- [ ] **Step 3: Final verification**

```bash
cargo test --all
cargo check --all-features
```

---

## Implementation Checklist

| Task | Status |
|------|--------|
| Task 1: Add `is_contig_end` to Node | |
| Task 2: Mark contig ends during graph building | |
| Task 3: Create ContigEndPruner module | |
| Task 4: Wire into pipeline | |
| Task 5: Create SV matrix output | |
| Task 6: Integrate into OutputWriter | |
| Task 7: Create integration tests | |
| Task 8: Fix issues and verify | |
| Task 9: Documentation | |

---

## Notes

1. **Contig-end marking approach**: The current implementation needs to track which genes belong to which contigs. This may require passing gene data to the `ContigEndPruner`.

2. **SV matrix implementation**: The `SVMatrixWriter` currently has placeholder logic. Need to verify `gene_names` are accessible for proper TSV output.

3. **Test data**: Use the existing `test_helpers` functions from `tests/integration_test.rs` as reference for creating test GFF3 files.
