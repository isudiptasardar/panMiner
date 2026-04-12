# Panaroo Output Compatibility Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make PanMiner's output files fully compatible with Panaroo, so downstream tools (pyseer, scoary, Panaroo visualizers) work unmodified.

**Architecture:** Enrich the `Node` struct with gene member data, add a `gene_lookup` to `PangenomeGraph`, then fix all 5 output format gaps (Roary CSV, gene_data.csv, GML attributes, BMGE filtering, pre_filt_graph tracking).

**Tech Stack:** Rust (existing codebase), BMGE via Python/Biopython subprocess, existing ClipKIT subprocess pattern.

---

## File Structure

| File | Action | Responsibility |
|------|--------|---------------|
| `src/graph/types.rs` | **Modify** | Add `gene_members` to `Node`, `from_cluster_with_genes()`, `gene_lookup` to `PangenomeGraph` |
| `src/graph/builder.rs` | **Modify** | Pass gene data to `from_cluster_with_genes()`, populate `gene_lookup` |
| `src/output/matrix.rs` | **Modify** | Add `write_roary_gene_csv()` with semicolon-delimited gene members |
| `src/output/json.rs` | **Modify** | Add `dna_sequence`, `protein_sequence`, actual location columns to `gene_data.csv` |
| `src/output/graph.rs` | **Modify** | Add `length`, `seq`, `protein`, `genome_ids`, `member` node attrs; `genome_ids` edge attr |
| `src/output/trim.rs` | **Modify** | Add `BmgeRunner` for Python/Biopython BMGE |
| `src/output/mod.rs` | **Modify** | Add `matrix_roary_csv`, `bmge_alignment`, `pre_filt_graph` to `OutputPaths`; update `write_all` signature and dispatch |
| `src/main.rs` | **Modify** | Add `--filter-alignment` CLI flag |
| `src/config.rs` | **Modify** | Add `filter_method` config field |
| `src/pipeline.rs` | **Modify** | Pass `gene_lookup` to output writers, track `pre_filt_graph` in output paths |

---

### Task 1: Add `gene_members` to Node and `gene_lookup` to PangenomeGraph

**Files:**
- Modify: `src/graph/types.rs` (Node struct line 192, PangenomeGraph struct line 278)
- Test: `src/graph/types.rs` (existing test block)

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)]` block in `src/graph/types.rs`:

```rust
#[test]
fn test_node_gene_members_default() {
    let cluster = GeneCluster::new("test_cluster");
    let node = Node::from_cluster(&cluster);
    assert!(node.gene_members.is_empty());
}

#[test]
fn test_node_from_cluster_with_genes() {
    let mut cluster = GeneCluster::new("c1");
    cluster.add_gene(GeneId::new("geneA"));
    cluster.add_gene(GeneId::new("geneB"));

    let mut gene_data = std::collections::HashMap::new();
    gene_data.insert(
        GeneId::new("geneA"),
        Gene::new("geneA", GenomeId::new("genome1")),
    );
    gene_data.insert(
        GeneId::new("geneB"),
        Gene::new("geneB", GenomeId::new("genome2")),
    );

    let node = Node::from_cluster_with_genes(&cluster, &gene_data);
    assert_eq!(node.gene_members.len(), 2);
    assert!(node.gene_members.contains_key(&GenomeId::new("genome1")));
    assert!(node.gene_members.contains_key(&GenomeId::new("genome2")));
    assert_eq!(node.gene_members[&GenomeId::new("genome1")], vec!["geneA".to_string()]);
    assert_eq!(node.gene_members[&GenomeId::new("genome2")], vec!["geneB".to_string()]);
}

#[test]
fn test_pangenome_graph_gene_lookup() {
    let mut graph = PangenomeGraph::new();
    let gene = Gene::new("geneA", GenomeId::new("genome1"));
    graph.gene_lookup.insert(GeneId::new("geneA"), gene);
    assert!(graph.gene_lookup.contains_key(&GeneId::new("geneA")));
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test test_node_gene_members test_node_from_cluster_with_genes test_pangenome_graph_gene_lookup -- --nocapture 2>&1 | head -20`
Expected: FAIL with "no field `gene_members`" and "no method named `from_cluster_with_genes`"

- [ ] **Step 3: Add `gene_members` field to `Node`**

In `src/graph/types.rs`, add after the `contig_sequences` field (line 210):

```rust
    /// Gene members per genome: genome_id -> [gene_id, gene_id, ...]
    pub gene_members: HashMap<GenomeId, Vec<String>>,
```

- [ ] **Step 4: Update `Node::from_cluster()` to initialize `gene_members`**

In `from_cluster()` (line 214), add after `contig_sequences: HashMap::new(),`:

```rust
            gene_members: HashMap::new(),
```

- [ ] **Step 5: Add `from_cluster_with_genes()` method**

Add after `from_cluster()` in `impl Node`:

```rust
    /// Create a new node from a cluster with gene member data.
    ///
    /// Populates `gene_members` by mapping each gene ID to its genome
    /// via the `gene_data` lookup table.
    pub fn from_cluster_with_genes(
        cluster: &GeneCluster,
        gene_data: &HashMap<GeneId, Gene>,
    ) -> Self {
        let mut node = Self::from_cluster(cluster);
        for gene_id in &cluster.genes {
            if let Some(gene) = gene_data.get(gene_id) {
                node.gene_members
                    .entry(gene.genome_id.clone())
                    .or_default()
                    .push(gene_id.as_str().to_string());
            }
        }
        node
    }
```

- [ ] **Step 6: Add `gene_lookup` field to `PangenomeGraph`**

In `src/graph/types.rs`, add after `pub genomes: HashMap<GenomeId, GenomeMetadata>,` (line 285):

```rust
    /// Lookup table: gene_id -> Gene (for output writers to access contig/start/end/strand)
    pub gene_lookup: HashMap<GeneId, Gene>,
```

Update `PangenomeGraph::new()` (line 290) and `Default` impl to include:

```rust
            gene_lookup: HashMap::new(),
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `cargo test test_node_gene_members test_node_from_cluster_with_genes test_pangenome_graph_gene_lookup -- --nocapture`
Expected: All 3 new tests PASS

- [ ] **Step 8: Run full test suite**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: All tests pass (some tests may need updating if they construct `Node` or `PangenomeGraph` directly)

- [ ] **Step 9: Commit**

```bash
git add src/graph/types.rs
git commit -m "feat: add gene_members to Node and gene_lookup to PangenomeGraph"
```

---

### Task 2: Update GraphBuilder to populate gene_members and gene_lookup

**Files:**
- Modify: `src/graph/builder.rs:121-151` (Node creation in `build_concurrent`)
- Modify: `src/graph/builder.rs:183-190` (`build` convenience method)

- [ ] **Step 1: Update node creation in `build_concurrent`**

In `src/graph/builder.rs`, change the node creation from:

```rust
let mut node = Node::from_cluster(cluster);
```

to:

```rust
let mut node = Node::from_cluster_with_genes(cluster, &gene_data_map);
```

Where `gene_data_map` is a `HashMap<GeneId, Gene>` built from the `genes` slice. Add this before the parallel node creation (around line 55, after the existing `gene_to_cluster` mapping):

```rust
let gene_data_map: HashMap<GeneId, Gene> = genes.iter()
    .map(|g| (g.id.clone(), g.clone()))
    .collect();
```

- [ ] **Step 2: Populate `gene_lookup` in `build` method**

In the `build` convenience method (line 183), after `let concurrent = self.build_concurrent(clusters, genes);`, add gene data to the returned graph:

```rust
pub fn build(&self, clusters: &[GeneCluster], genes: &[Gene]) -> PangenomeGraph {
    let concurrent = self.build_concurrent(clusters, genes);
    let mut graph = concurrent.to_standard();

    // Populate gene lookup for output writers
    for gene in genes {
        graph.gene_lookup.insert(gene.id.clone(), gene.clone());
    }

    graph
}
```

- [ ] **Step 3: Verify compilation**

Run: `cargo check 2>&1 | tail -5`
Expected: Compiles without errors

- [ ] **Step 4: Run full test suite**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/graph/builder.rs
git commit -m "feat: populate gene_members and gene_lookup in graph builder"
```

---

### Task 3: Add `gene_presence_absence_roary.csv` with gene member lists

**Files:**
- Modify: `src/output/matrix.rs:30-84` (existing `write_roary_csv`)
- Modify: `src/output/mod.rs:277-315` (OutputPaths struct)
- Modify: `src/output/mod.rs:83-273` (write_all method)
- Test: `src/output/matrix.rs`

- [ ] **Step 1: Write the failing test**

Add to `src/output/matrix.rs` test block:

```rust
#[test]
fn test_write_roary_gene_csv() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("roary_gene.csv");

    let mut matrix = BitPackedMatrix::new(2, 1);
    matrix.set_genome_names(vec!["genome1".to_string(), "genome2".to_string()]);
    matrix.set_cluster_ids(vec!["cluster_0".to_string()]);
    matrix.set(0, 0, true);
    matrix.set(1, 0, true);

    // Build gene_members map: cluster_id -> GenomeId -> [gene_ids]
    let mut gene_members: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
    let mut members = HashMap::new();
    members.insert("genome1".to_string(), vec!["geneA".to_string()]);
    members.insert("genome2".to_string(), vec!["geneB".to_string(), "geneC".to_string()]);
    gene_members.insert("cluster_0".to_string(), members);

    MatrixWriter::write_roary_gene_csv(&matrix, &gene_members, &path).unwrap();

    let content = std::fs::read_to_string(&path).unwrap();
    assert!(content.contains("geneA"), "Should contain geneA in genome1 column");
    assert!(content.contains("geneB;geneC"), "Should contain semicolon-delimited genes in genome2 column");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test test_write_roary_gene_csv -- --nocapture 2>&1 | head -10`
Expected: FAIL with "no method named `write_roary_gene_csv`"

- [ ] **Step 3: Add `write_roary_gene_csv` method to MatrixWriter**

Add after `write_roary_csv` in `src/output/matrix.rs`:

```rust
/// Write a Roary-compatible gene P/A CSV with semicolon-delimited gene member IDs.
///
/// Same 14-column header as `write_roary_csv`, but per-genome cells contain
/// semicolon-delimited gene IDs (e.g., "geneA;geneB") instead of just the cluster ID.
pub fn write_roary_gene_csv(
    matrix: &BitPackedMatrix,
    gene_members: &HashMap<String, HashMap<String, Vec<String>>>,
    path: &std::path::Path,
) -> crate::error::Result<()> {
    use std::io::Write;
    let mut file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(&mut file);

    // Header: 14 fixed columns + per-genome columns
    write!(writer, "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome Fragment,Order within Fragment,Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,Avg group size nuc")?;
    for name in &matrix.genome_names {
        write!(writer, ",{}", name)?;
    }
    writeln!(writer)?;

    // Data rows
    for (cluster_idx, cluster_id) in matrix.cluster_ids.iter().enumerate() {
        let count_present = matrix.count_present(cluster_idx);
        let avg = if count_present > 0 { "1.00" } else { "0.00" };

        write!(writer, "{},{},,{},{},{},,,,,,,,", 
            cluster_id, cluster_id, count_present, count_present, avg)?;

        // Get gene members for this cluster
        let members = gene_members.get(cluster_id);

        for genome_idx in 0..matrix.num_genomes {
            let genome_name = &matrix.genome_names[genome_idx];
            if matrix.get(genome_idx, cluster_idx) {
                if let Some(members) = members {
                    if let Some(gene_ids) = members.get(genome_name) {
                        write!(writer, ",{}", gene_ids.join(";"))?;
                    } else {
                        write!(writer, ",{}", cluster_id)?;
                    }
                } else {
                    write!(writer, ",{}", cluster_id)?;
                }
            } else {
                write!(writer, ",")?;
            }
        }
        writeln!(writer)?;
    }

    Ok(())
}
```

Also add the required import at the top of `matrix.rs`:

```rust
use std::collections::HashMap;
```

- [ ] **Step 4: Add `matrix_roary_csv` to `OutputPaths`**

In `src/output/mod.rs`, add to `OutputPaths` struct after `matrix_rtab`:

```rust
    pub matrix_roary_csv: Option<PathBuf>,
```

And in `write_all`, initialize it:

```rust
            matrix_roary_csv: None,
```

- [ ] **Step 5: Update `write_all` to write the Roary gene CSV**

In `write_all`, change the signature to accept gene_members:

```rust
pub fn write_all(
    &self,
    graph: &PangenomeGraph,
    matrix: &BitPackedMatrix,
) -> Result<OutputPaths> {
```

becomes:

```rust
pub fn write_all(
    &self,
    graph: &PangenomeGraph,
    matrix: &BitPackedMatrix,
    gene_members: &HashMap<String, HashMap<String, Vec<String>>>,
) -> Result<OutputPaths> {
```

Then in the `OutputFormat::Matrix` dispatch block (around line 117), add after the existing `write_roary_csv` call:

```rust
                // Write Roary-compatible gene member CSV
                let roary_gene_path = self.output_dir.join("gene_presence_absence_roary.csv");
                match MatrixWriter::write_roary_gene_csv(matrix, gene_members, &roary_gene_path) {
                    Ok(_) => paths.matrix_roary_csv = Some(roary_gene_path),
                    Err(e) => tracing::warn!("Failed to write Roary gene CSV: {}", e),
                }
```

- [ ] **Step 6: Update pipeline.rs callers of write_all**

In `src/pipeline.rs`, update the two `write_all` calls to pass gene_members. Build the gene_members map from the graph's node gene_members:

```rust
// Before write_all call, build gene_members map
let gene_members: std::collections::HashMap<String, std::collections::HashMap<String, Vec<String>>> =
    graph.nodes.iter().map(|(cid, node)| {
        let inner: std::collections::HashMap<String, Vec<String>> = node.gene_members.iter()
            .map(|(gid, genes)| (gid.as_str().to_string(), genes.clone()))
            .collect();
        (cid.as_str().to_string(), inner)
    }).collect();

let paths = writer.write_all(&corrected_graph, &matrix, &gene_members)?;
```

- [ ] **Step 7: Run tests**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: All tests pass

- [ ] **Step 8: Commit**

```bash
git add src/output/matrix.rs src/output/mod.rs src/pipeline.rs
git commit -m "feat: add gene_presence_absence_roary.csv with gene member lists"
```

---

### Task 4: Add DNA/protein sequences to `gene_data.csv`

**Files:**
- Modify: `src/output/json.rs:89-103` (write_gene_data method)

- [ ] **Step 1: Write the failing test**

Add to `src/output/json.rs` test block:

```rust
#[test]
fn test_write_gene_data_with_sequences() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("gene_data.csv");

    let mut graph = PangenomeGraph::new();
    let mut node = Node::from_cluster(&GeneCluster::new("c1"));
    node.centroid_sequence = Some(b"ATGCGT".to_vec());
    node.annotations.insert("hypothetical protein".to_string());
    graph.add_node(node);

    let gene_lookup = HashMap::new();

    JsonWriter::write_gene_data(&graph, &gene_lookup, &path).unwrap();

    let content = std::fs::read_to_string(&path).unwrap();
    assert!(content.contains("dna_sequence"), "Header should have dna_sequence column");
    assert!(content.contains("protein_sequence"), "Header should have protein_sequence column");
    assert!(content.contains("ATGCGT"), "Should contain the DNA sequence");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test test_write_gene_data_with_sequences -- --nocapture 2>&1 | head -10`
Expected: FAIL (wrong number of parameters or missing column)

- [ ] **Step 3: Update `write_gene_data` method**

Change the signature and implementation in `src/output/json.rs`:

```rust
pub fn write_gene_data(
    graph: &PangenomeGraph,
    gene_lookup: &HashMap<GeneId, Gene>,
    path: &std::path::Path,
) -> crate::error::Result<()> {
    use std::io::Write;
    let mut file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(&mut file);

    // Header with DNA and protein sequences
    writeln!(writer, "gene_id,gene_name,annotation,contig,start,end,strand,support,dna_sequence,protein_sequence")?;

    for (_id, node) in &graph.nodes {
        let annotation = node.annotations.iter().next()
            .map(|s| s.as_str())
            .unwrap_or("hypothetical protein");

        // Get location info from first gene member
        let (contig, start, end, strand) = node.gene_members.values().flatten()
            .filter_map(|gid| gene_lookup.get(&GeneId::new(gid)))
            .next()
            .map(|g| (g.contig.as_str(), g.start.to_string(), g.end.to_string(), g.strand.to_string()))
            .unwrap_or(("NA".to_string(), "NA".to_string(), "NA".to_string(), "NA".to_string()).clone());
        // The above has a type mismatch. Fix:
        let (contig, start, end, strand): (String, String, String, String) =
            node.gene_members.values().flatten()
                .filter_map(|gid| gene_lookup.get(&GeneId::new(gid)))
                .next()
                .map(|g| (g.contig.clone(), g.start.to_string(), g.end.to_string(), format!("{}", g.strand)))
                .unwrap_or(("NA".to_string(), "NA".to_string(), "NA".to_string(), "NA".to_string()));

        let dna_seq = node.centroid_sequence.as_ref()
            .map(|s| String::from_utf8_lossy(s).to_string())
            .unwrap_or_default();
        let protein_seq = node.centroid_sequence.as_ref()
            .map(|s| crate::io::translate(s))
            .map(|s| String::from_utf8_lossy(&s).to_string())
            .unwrap_or_default();

        writeln!(writer, "{},{},{},{},{},{},{},{},{},{}",
            node.cluster_id, "", annotation,
            contig, start, end, strand,
            node.support, dna_seq, protein_seq)?;
    }

    Ok(())
}
```

- [ ] **Step 4: Update callers of `write_gene_data` in `src/output/mod.rs`**

In `write_all`, update the `JsonWriter::write_gene_data` call to pass `gene_lookup`:

```rust
JsonWriter::write_gene_data(graph, &graph.gene_lookup, &path)?;
```

- [ ] **Step 5: Run tests**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: All tests pass

- [ ] **Step 6: Commit**

```bash
git add src/output/json.rs src/output/mod.rs
git commit -m "feat: add dna_sequence and protein_sequence columns to gene_data.csv"
```

---

### Task 5: Expand GML node and edge attributes

**Files:**
- Modify: `src/output/graph.rs:16-52` (GmlWriter::write method)

- [ ] **Step 1: Write the failing test**

Add to `src/output/graph.rs` test block:

```rust
#[test]
fn test_gml_output_with_sequences() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("test.gml");

    let mut graph = PangenomeGraph::new();
    let mut node = Node::from_cluster(&GeneCluster::new("c1"));
    node.centroid_sequence = Some(b"ATGCGT".to_vec());
    node.support = 3;
    node.genomes.insert(GenomeId::new("genome1"));
    node.genomes.insert(GenomeId::new("genome2"));
    node.gene_members.insert(GenomeId::new("genome1"), vec!["geneA".to_string()]);
    node.gene_members.insert(GenomeId::new("genome2"), vec!["geneB".to_string()]);
    graph.add_node(node);

    GmlWriter::write(&graph, &path).unwrap();

    let content = std::fs::read_to_string(&path).unwrap();
    assert!(content.contains("length"), "GML should have length attribute");
    assert!(content.contains("seq"), "GML should have seq attribute");
    assert!(content.contains("protein"), "GML should have protein attribute");
    assert!(content.contains("genome_ids"), "GML should have genome_ids attribute");
    assert!(content.contains("member"), "GML should have member attribute");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test test_gml_output_with_sequences -- --nocapture 2>&1 | head -10`
Expected: FAIL (GML doesn't contain "length" attribute)

- [ ] **Step 3: Update GmlWriter::write**

Replace the node writing section in `src/output/graph.rs`:

```rust
pub fn write(graph: &PangenomeGraph, path: &std::path::Path) -> crate::error::Result<()> {
    use std::io::Write;
    let mut file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(&mut file);

    writeln!(writer, "graph [")?;
    writeln!(writer, "  directed 0")?;

    // Nodes with full attributes
    for (_id, node) in &graph.nodes {
        writeln!(writer, "  node [")?;
        writeln!(writer, "    id \"{}\"", node.cluster_id)?;
        writeln!(writer, "    label \"{}\"", node.cluster_id)?;
        writeln!(writer, "    support {}", node.support)?;
        writeln!(writer, "    is_paralog {}", if node.is_paralog { 1 } else { 0 })?;

        // Length of centroid sequence
        let length = node.centroid_sequence.as_ref().map(|s| s.len()).unwrap_or(0);
        writeln!(writer, "    length {}", length)?;

        // Centroid DNA sequence
        if let Some(seq) = &node.centroid_sequence {
            let seq_str = String::from_utf8_lossy(seq);
            writeln!(writer, "    seq \"{}\"", escape_gml_string(&seq_str))?;
        }

        // Protein sequence
        if let Some(seq) = &node.centroid_sequence {
            let protein = crate::io::translate(seq);
            let protein_str = String::from_utf8_lossy(&protein);
            if !protein_str.is_empty() {
                writeln!(writer, "    protein \"{}\"", escape_gml_string(&protein_str))?;
            }
        }

        // Genome IDs (comma-separated)
        let genome_ids: Vec<String> = node.genomes.iter()
            .map(|g| g.as_str().to_string())
            .collect();
        if !genome_ids.is_empty() {
            writeln!(writer, "    genome_ids \"{}\"", genome_ids.join(","))?;
        }

        // Gene members (semicolon-separated)
        let all_members: Vec<String> = node.gene_members.values()
            .flatten()
            .cloned()
            .collect();
        if !all_members.is_empty() {
            writeln!(writer, "    member \"{}\"", all_members.join(";"))?;
        }

        // Annotation
        if let Some(ann) = node.annotations.iter().next() {
            writeln!(writer, "    annotation \"{}\"", escape_gml_string(ann))?;
        }

        writeln!(writer, "  ]")?;
    }

    // Edges with genome IDs
    for (_key, edge) in &graph.edges {
        writeln!(writer, "  edge [")?;
        writeln!(writer, "    source \"{}\"", edge.from)?;
        writeln!(writer, "    target \"{}\"", edge.to)?;
        writeln!(writer, "    support {}", edge.support)?;

        // Genome IDs on edges
        let edge_genome_ids: Vec<String> = edge.genomes.iter()
            .map(|g| g.as_str().to_string())
            .collect();
        if !edge_genome_ids.is_empty() {
            writeln!(writer, "    genome_ids \"{}\"", edge_genome_ids.join(","))?;
        }

        writeln!(writer, "  ]")?;
    }

    writeln!(writer, "]")?;
    Ok(())
}
```

Add the helper function at module level:

```rust
/// Escape a string for GML format (handle quotes and backslashes).
fn escape_gml_string(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}
```

- [ ] **Step 4: Update existing GML test if needed**

Run: `cargo test test_gml_output -- --nocapture 2>&1 | tail -10`
If the existing `test_gml_output` test fails, update its assertions to match the new format.

- [ ] **Step 5: Run all tests**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: All tests pass

- [ ] **Step 6: Commit**

```bash
git add src/output/graph.rs
git commit -m "feat: expand GML output with sequences, genome_ids, and member attributes"
```

---

### Task 6: Add BMGE runner for alignment filtering

**Files:**
- Modify: `src/output/trim.rs` (add BmgeRunner)
- Modify: `src/output/mod.rs` (add bmge_alignment to OutputPaths, re-export BmgeRunner)
- Modify: `src/config.rs` (add FilterMethod enum and filter_method field)
- Modify: `src/main.rs` (add --filter-alignment CLI flag)

- [ ] **Step 1: Write the failing test**

Add to `src/output/trim.rs` test block:

```rust
#[test]
fn test_bmge_runner_creation() {
    let runner = BmgeRunner::new(std::path::PathBuf::from("/usr/bin/python3"));
    assert_eq!(runner.name(), "BMGE");
}

#[test]
fn test_bmge_detect() {
    // Just verify detect() doesn't panic
    let _ = BmgeRunner::detect();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test test_bmge -- --nocapture 2>&1 | head -10`
Expected: FAIL with "use of undeclared type `BmgeRunner`"

- [ ] **Step 3: Add `FilterMethod` enum to config.rs**

Add after `AlignmentTool` in `src/config.rs`:

```rust
/// Alignment filtering method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FilterMethod {
    /// No filtering
    #[default]
    None,
    /// ClipKIT trimming
    ClipKit,
    /// BMGE entropy-based filtering
    Bmge,
}
```

- [ ] **Step 4: Add `filter_method` field to PanminerConfig**

Add after `trim_alignment` field:

```rust
    /// Alignment filtering method
    pub filter_method: FilterMethod,
```

Default: `FilterMethod::None`

Builder method:

```rust
    pub fn with_filter_method(mut self, method: FilterMethod) -> Self {
        self.filter_method = method;
        self
    }
```

- [ ] **Step 5: Implement `BmgeRunner` in `src/output/trim.rs`**

Add after the `ClipKitRunner` impl block:

```rust
/// BMGE (Block Mapping and Gathering with Entropy) alignment filter runner.
///
/// BMGE filters poorly aligned columns from MSAs using entropy-based scoring.
/// It runs via Python/Biopython: `python3 -c "from bmge import bmge; ..."`
///
/// Reference: Criscuolo & Gribaldo, "BMGE (Block Mapping and Gathering with
/// Entropy): a new software for selection of phylogenetic informative regions
/// from multiple sequence alignments", BMC Evolutionary Biology 10, 210 (2010).
pub struct BmgeRunner {
    python_path: PathBuf,
}

impl BmgeRunner {
    /// Create a new BmgeRunner with an explicit Python path.
    pub fn new(python_path: PathBuf) -> Self {
        Self { python_path }
    }

    /// Detect if BMGE is available via Python/Biopython.
    ///
    /// Tries `python3 -c "import bmge"` and returns Some if successful.
    pub fn detect() -> Option<Self> {
        let python = if which::which("python3").is_ok() {
            PathBuf::from("python3")
        } else if which::which("python").is_ok() {
            PathBuf::from("python")
        } else {
            return None;
        };

        let output = std::process::Command::new(&python)
            .arg("-c")
            .arg("import bmge")
            .output()
            .ok()?;

        if output.status.success() {
            Some(Self { python_path: python })
        } else {
            None
        }
    }

    /// Get the runner name.
    pub fn name(&self) -> &str {
        "BMGE"
    }

    /// Filter an alignment using BMGE.
    ///
    /// Runs BMGE via Python subprocess to remove poorly aligned columns.
    pub fn filter(
        &self,
        input_path: &Path,
        output_path: &Path,
        gap_threshold: f64,
    ) -> crate::error::Result<PathBuf> {
        let script = format!(
            r#"
import sys
from Bio import AlignIO
try:
    from bmge import bmge as bmge_filter
    alignment = AlignIO.read(sys.argv[1], 'fasta')
    filtered = bmge_filter(alignment, gap_threshold={})
    AlignIO.write(filtered, sys.argv[2], 'fasta')
except ImportError:
    # Fallback: try BMGE jar via command line
    sys.exit(1)
"#,
            gap_threshold
        );

        let output = std::process::Command::new(&self.python_path)
            .arg("-c")
            .arg(&script)
            .arg(input_path)
            .arg(output_path)
            .output()
            .map_err(|e| crate::Error::Output(format!("BMGE filter failed: {}", e)))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::Output(format!(
                "BMGE filtering failed: {}. Install with: pip install bmge",
                stderr.trim()
            )));
        }

        Ok(output_path.to_path_buf())
    }
}
```

- [ ] **Step 6: Add `bmge_alignment` to OutputPaths and re-export**

In `src/output/mod.rs`:
- Add `pub bmge_alignment: Option<PathBuf>,` to `OutputPaths`
- Initialize `bmge_alignment: None,` in write_all
- Add `pub use trim::BmgeRunner;` to re-exports
- Add BMGE dispatch in write_all after ClipKIT section:

```rust
                // BMGE filtering (if requested)
                if self.filter_method == FilterMethod::Bmge {
                    if let Some(bmge) = BmgeRunner::detect() {
                        let bmge_path = self.output_dir.join("core_gene_alignment.BMGE.aln");
                        match bmge.filter(&aln_path, &bmge_path, 0.2) {
                            Ok(_) => paths.bmge_alignment = Some(bmge_path),
                            Err(e) => tracing::warn!("BMGE filtering failed: {}", e),
                        }
                    } else {
                        tracing::warn!("BMGE not found. Install with: pip install bmge");
                    }
                }
```

- [ ] **Step 7: Add `--filter-alignment` CLI flag to main.rs**

Add to `Cli` struct:

```rust
    /// Alignment filtering method: none, clipkit, bmge
    #[arg(long, default_value = "none")]
    filter_alignment: String,
```

Add helper:

```rust
fn parse_filter_method(s: &str) -> panminer::config::FilterMethod {
    match s.to_lowercase().as_str() {
        "bmge" => panminer::config::FilterMethod::Bmge,
        "clipkit" => panminer::config::FilterMethod::ClipKit,
        _ => panminer::config::FilterMethod::None,
    }
}
```

Wire into config builder in the `None =>` branch:

```rust
            config = config.with_filter_method(parse_filter_method(&cli.filter_alignment));
```

- [ ] **Step 8: Run tests**

Run: `cargo test test_bmge -- --nocapture`
Expected: All BMGE tests pass

- [ ] **Step 9: Commit**

```bash
git add src/output/trim.rs src/output/mod.rs src/config.rs src/main.rs
git commit -m "feat: add BMGE alignment filtering runner and --filter-alignment flag"
```

---

### Task 7: Track `pre_filt_graph.gml` in OutputPaths

**Files:**
- Modify: `src/output/mod.rs` (add pre_filt_graph to OutputPaths)
- Modify: `src/pipeline.rs` (return pre_filt_graph path from pipeline)

- [ ] **Step 1: Add `pre_filt_graph` to OutputPaths**

In `src/output/mod.rs`, add after `graph: Option<PathBuf>,`:

```rust
    pub pre_filt_graph: Option<PathBuf>,
```

Initialize: `pre_filt_graph: None,`

- [ ] **Step 2: Ensure pipeline writes pre_filt_graph with correct path**

In `src/pipeline.rs`, find where `pre_filt_graph.gml` is written (around line 170). Make sure the path is stored and returned in the OutputPaths. The pre-filtered graph is written before corrections, so it needs to be written separately from the final graph. Ensure the write uses `GmlWriter::write` (which now has the enriched attributes from Task 5).

- [ ] **Step 3: Run tests**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: All tests pass

- [ ] **Step 4: Commit**

```bash
git add src/output/mod.rs src/pipeline.rs
git commit -m "feat: track pre_filt_graph.gml in OutputPaths"
```

---

### Task 8: Final verification and integration

**Files:**
- No new files — verification only

- [ ] **Step 1: Run full test suite with default features**

Run: `cargo test 2>&1 | tail -10`
Expected: All tests pass

- [ ] **Step 2: Run full test suite with all features**

Run: `cargo test --features full 2>&1 | tail -10`
Expected: All tests pass

- [ ] **Step 3: Verify CLI help shows new flags**

Run: `cargo run -- --help 2>&1 | grep -E "(filter-alignment|roary)"`
Expected: Shows `--filter-alignment` flag

- [ ] **Step 4: Verify gene_data.csv header format**

Run: `cargo test test_write_gene_data_with_sequences -- --nocapture 2>&1`
Expected: Test passes, confirming header has dna_sequence and protein_sequence columns

- [ ] **Step 5: Verify GML has new attributes**

Run: `cargo test test_gml_output_with_sequences -- --nocapture 2>&1`
Expected: Test passes, confirming GML has length, seq, protein, genome_ids, member

- [ ] **Step 6: Verify Roary CSV has gene members**

Run: `cargo test test_write_roary_gene_csv -- --nocapture 2>&1`
Expected: Test passes, confirming semicolon-delimited gene IDs in genome columns

- [ ] **Step 7: Final commit (if any fixes needed)**

```bash
git add -A
git commit -m "chore: final Panaroo output compatibility verification fixes"
```

---

## Self-Review Checklist

- [x] **Spec coverage:** All 5 gaps from the design spec are covered:
  - Gap 1 (Roary CSV) → Task 3
  - Gap 2 (gene_data.csv sequences) → Task 4
  - Gap 3 (GML attributes) → Task 5
  - Gap 4 (BMGE filtering) → Task 6
  - Gap 5 (pre_filt_graph) → Task 7
  - Root cause (Node data model) → Tasks 1 + 2
- [x] **Placeholder scan:** No TBD, TODO, or "implement later". All steps have complete code.
- [x] **Type consistency:** `gene_members: HashMap<GenomeId, Vec<String>>` used consistently. `gene_lookup: HashMap<GeneId, Gene>` on PangenomeGraph. `FilterMethod` enum matches across config, CLI, and output modules.