# Panaroo Feature Parity — Phase 2: New Subcommands

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `integrate`, `extract-gene`, `msa`, and GFF output subcommands to match Panaroo's CLI coverage.

**Architecture:** All four are new code with no existing module changes. Integrate is the most complex (loads GML, parses GFF, clusters against centroids, adds genes to graph, runs corrections). Extract-gene and MSA are simple lookups/wrappers. GFF output reconstructs GFF3 from the corrected graph.

**Tech Stack:** Rust, clap, existing pipeline modules

---

### Task 1: Add `extract-gene` subcommand

**Files:**
- Create: `src/io/extract_gene.rs`
- Modify: `src/main.rs` (add ExtractGene subcommand)
- Modify: `src/lib.rs` (add module)

- [ ] **Step 1: Create the extract_gene module**

Create `src/io/extract_gene.rs`:

```rust
//! Extract individual gene sequences from a pangenome output directory.

use std::path::Path;
use std::collections::HashMap;
use crate::error::Result;

/// Extract sequences for a given cluster from PanMiner output.
pub fn extract_gene(
    output_dir: &Path,
    cluster_id: &str,
    output_path: &Path,
    protein: bool,
) -> Result<()> {
    // 1. Load final_graph.gml to find the node by cluster_id
    // 2. Load gene_data.csv to get per-gene sequences
    // 3. Filter to genes belonging to the cluster
    // 4. Write FASTA to output_path
    //    If protein=true, write protein sequences; otherwise DNA
    todo!("implement after GML reader is available")
}
```

- [ ] **Step 2: Add ExtractGene subcommand to CLI**

In `src/main.rs`, add to the `Commands` enum:

```rust
/// Extract gene sequences for a specific cluster
#[command(name = "extract-gene")]
ExtractGene {
    /// PanMiner output directory
    #[arg(short = 'i', long)]
    input: PathBuf,

    /// Cluster ID to extract
    #[arg(long)]
    cluster: String,

    /// Output FASTA file
    #[arg(short = 'o', long, default_value = "extracted_genes.fasta")]
    output: PathBuf,

    /// Extract protein sequences instead of DNA
    #[arg(long)]
    protein: bool,
},
```

Wire in the match arm:

```rust
Commands::ExtractGene { input, cluster, output, protein } => {
    panminer::io::extract_gene::extract_gene(&input, &cluster, &output, protein)?;
}
```

- [ ] **Step 3: Add module to lib.rs**

```rust
pub mod io {
    // ... existing modules ...
    pub mod extract_gene;
}
```

- [ ] **Step 4: Implement extract_gene using existing gene_data.csv**

Read `gene_data.csv` (which has gene_id, genome_id, cluster_id, dna_sequence, protein_sequence columns). Filter rows matching the cluster ID. Write matching sequences as FASTA.

- [ ] **Step 5: Write test**

```rust
#[test]
fn test_extract_gene_finds_cluster() {
    // Create a small gene_data.csv with known clusters
    // Call extract_gene with a cluster ID
    // Verify the output FASTA contains the correct sequences
}
```

- [ ] **Step 6: Run tests and commit**

```bash
cargo test --features full
git add src/io/extract_gene.rs src/main.rs src/lib.rs
git commit -m "feat: add extract-gene subcommand for cluster sequence extraction"
```

---

### Task 2: Add GFF output writer

**Files:**
- Create: `src/output/gff.rs`
- Modify: `src/output/mod.rs` (register module)
- Modify: `src/config.rs` (add Gff to OutputFormat enum)
- Modify: `src/main.rs` (add `gff` to formats argument)

- [ ] **Step 1: Create the GFF writer module**

Create `src/output/gff.rs`:

```rust
//! Generate GFF3 output files from the corrected pangenome graph.

use std::path::Path;
use std::collections::HashMap;
use crate::error::Result;
use crate::graph::PangenomeGraph;

/// Write GFF3 files for each genome from the pangenome graph.
pub fn write_gff_files(
    graph: &PangenomeGraph,
    output_dir: &Path,
) -> Result<Vec<std::path::PathBuf>> {
    // 1. Group nodes by genome_id (from gene_members)
    // 2. For each genome:
    //    a. Sort nodes by contig and position (from gene_data)
    //    b. Write GFF3 header: ##gff-version 3
    //    c. Write CDS features with:
    //       - seqid = contig name
    //       - source = "panminer"
    //       - type = "CDS"
    //       - start, end = gene coordinates
    //       - score = "."
    //       - strand = +/-
    //       - phase = "0"
    //       - attributes: ID=gene_id;cluster=cluster_id;annotation=...
    //    d. Write ##FASTA section with contig sequences
    // 3. Return list of written file paths
    todo!("implement after graph iteration is available")
}
```

- [ ] **Step 2: Add `Gff` variant to OutputFormat enum**

In `src/config.rs`:

```rust
pub enum OutputFormat {
    // ... existing variants ...
    Gff,
}
```

Update `parse_formats` to include `"gff"`.

- [ ] **Step 3: Register in output module**

In `src/output/mod.rs`:

```rust
pub mod gff;
// Add to OutputWriter::write_all:
if formats.contains(&OutputFormat::Gff) {
    gff::write_gff_files(&graph, &paths.output_dir)?;
}
```

- [ ] **Step 4: Implement GFF3 writing**

Use the `gene_members` HashMap on each Node to reconstruct per-genome GFF records. Each gene entry gets:
- `seqid`: contig name from gene data
- `source`: "panminer"
- `type`: "CDS"
- `start`, `end`: from gene data
- `strand`: from gene data
- `attributes`: `ID=<gene_id>;cluster=<cluster_id>`

- [ ] **Step 5: Write test and commit**

```bash
cargo test --features full
git add src/output/gff.rs src/output/mod.rs src/config.rs src/main.rs
git commit -m "feat: add GFF3 output writer for corrected pangenome annotations"
```

---

### Task 3: Add `msa` subcommand for standalone alignment

**Files:**
- Modify: `src/main.rs` (add Msa subcommand)
- Modify: `src/clustering/alignment_traits.rs` (add standalone entry point)

- [ ] **Step 1: Add Msa subcommand to CLI**

In `src/main.rs`, add to the `Commands` enum:

```rust
/// Run multiple sequence alignment on pangenome output
#[command(name = "msa")]
Msa {
    /// PanMiner output directory containing final_graph.gml
    #[arg(short = 'i', long)]
    input: PathBuf,

    /// Output directory for alignment files
    #[arg(short = 'o', long, default_value = "msa_output")]
    output: PathBuf,

    /// Alignment mode: core or pan
    #[arg(long, default_value = "core")]
    mode: String,

    /// Alignment tool: mafft, clustal, prank
    #[arg(long, default_value = "mafft")]
    aligner: String,

    /// Number of threads
    #[arg(short = 't', long, default_value = "1")]
    threads: usize,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
},
```

- [ ] **Step 2: Implement Msa handler**

The handler loads `final_graph.gml` and `gene_data.csv`, identifies core/accessory genes based on the `--mode` argument, and runs alignment using existing `AlignmentRunner` implementations.

```rust
Commands::Msa { input, output, mode, aligner, threads, verbose } => {
    // 1. Load graph from GML
    // 2. Load gene_data.csv for sequences
    // 3. Select genes based on mode (core = support >= 99% of genomes)
    // 4. Create aligner based on --aligner flag
    // 5. Run alignment per gene
    // 6. Write alignment files to output directory
}
```

- [ ] **Step 3: Write test and commit**

```bash
cargo test --features full
git add src/main.rs
git commit -m "feat: add msa subcommand for standalone multiple sequence alignment"
```

---

### Task 4: Add `integrate` subcommand

**Files:**
- Create: `src/io/integrate.rs`
- Modify: `src/main.rs` (add Integrate subcommand)
- Modify: `src/lib.rs` (add module)

- [ ] **Step 1: Create the integrate module**

Create `src/io/integrate.rs`:

```rust
//! Integrate a new genome into an existing pangenome graph.

use crate::error::Result;
use crate::graph::PangenomeGraph;
use crate::config::PanminerConfig;
use std::path::Path;

/// Add a single GFF file to an existing PanMiner pangenome.
pub fn integrate_genome(
    existing_dir: &Path,
    new_gff: &Path,
    output_dir: &Path,
    config: &PanminerConfig,
) -> Result<()> {
    // 1. Load final_graph.gml from existing_dir
    // 2. Parse new GFF file
    // 3. Cluster new genes against existing centroids using MMseqs2 easy-search
    //    or CPU fallback with length filtering
    // 4. Add matching genes to existing nodes (increment support, add gene_members)
    // 5. Create new nodes for genes below identity threshold
    // 6. Build adjacency edges for new genes on the same contig
    // 7. Run correction passes on the modified graph
    // 8. Write updated output files
    todo!("implement")
}
```

- [ ] **Step 2: Add Integrate subcommand to CLI**

```rust
/// Integrate a new genome into an existing pangenome
#[command(name = "integrate")]
Integrate {
    /// Existing PanMiner output directory
    #[arg(long)]
    graph: PathBuf,

    /// New GFF3 file to integrate
    #[arg(long)]
    input: PathBuf,

    /// Output directory for updated pangenome
    #[arg(short = 'o', long, default_value = "integrated_output")]
    output: PathBuf,

    /// Identity threshold for matching new genes (0.5-1.0)
    #[arg(long, default_value = "0.98")]
    identity: f32,

    /// Number of threads (0 = auto-detect)
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
},
```

- [ ] **Step 3: Implement core integration logic**

The key algorithm:
1. Load existing GML graph into `PangenomeGraph`
2. Parse new GFF with `GffParser`
3. For each new gene, compute identity against existing centroids using `CpuClusterer::sequence_identity` or MMseqs2 search
4. If match above threshold: add gene to existing node's `gene_members`, increment `support`
5. If below threshold: create new node with new cluster
6. Build edges between consecutive genes on the same contig in the new genome
7. Run corrections (paralog resolution, contamination removal, etc.)

- [ ] **Step 4: Write test and commit**

```bash
cargo test --features full
git add src/io/integrate.rs src/main.rs src/lib.rs
git commit -m "feat: add integrate subcommand for incremental genome addition"
```