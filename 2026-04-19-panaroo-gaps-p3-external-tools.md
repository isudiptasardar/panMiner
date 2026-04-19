# Panaroo Feature Parity — Phase 3: External Tool Wiring

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire Pangrowth (pangenome openness), add abundance visualization (HTML/D3.js), and add Orphos/Prodigal gene calling (Rust crate).

**Architecture:** Pangrowth and Orphos are external tools wired as subprocesses/crates. Abundance visualization reuses the D3.js HTML pattern from `qc_viz.rs`. All three are independent features.

**Tech Stack:** Rust, orphos-core crate (feature-gated), subprocess for pangrowth, D3.js for visualization

---

### Task 1: Add Pangrowth pangenome openness subprocess

**Files:**
- Create: `src/downstream/evolution/pangrowth.rs`
- Modify: `src/downstream/evolution/mod.rs` (register module)
- Modify: `src/main.rs` (add `--pangrowth` flag to analyze subcommand)
- Modify: `src/downstream/mod.rs` (add re-export)

- [ ] **Step 1: Create PangrowthRunner module**

Create `src/downstream/evolution/pangrowth.rs`:

```rust
//! Pangrowth subprocess runner for pangenome openness estimation.
//!
//! Computes exact pangenome growth/core curves and fits Heaps' law alpha
//! to classify pangenomes as open or closed.
//! Reference: Parmigiani, Wittler, Stoye (2024) PCI Comp Biol.

use std::path::{Path, PathBuf};
use std::process::Command;
use crate::error::{Error, Result};

/// Pangenome openness classification.
#[derive(Debug, Clone, PartialEq)]
pub enum OpennessClassification {
    Open,    // alpha > 0: new genes continue indefinitely
    Closed,  // alpha <= 0: pangenome size converges
}

/// Results from Pangrowth analysis.
#[derive(Debug, Clone)]
pub struct PangrowthResult {
    /// Heaps' law alpha parameter (exponent)
    pub alpha: f64,
    /// Heaps' law kappa parameter (coefficient)
    pub kappa: f64,
    /// Openness classification
    pub classification: OpennessClassification,
    /// Growth curve data: (n_genomes, expected_pangenome_size)
    pub growth_curve: Vec<(usize, f64)>,
    /// Core curve data: (n_genomes, expected_core_size)
    pub core_curve: Vec<(usize, f64)>,
}

pub struct PangrowthRunner {
    path: PathBuf,
}

impl PangrowthRunner {
    /// Detect Pangrowth on the system PATH.
    pub fn detect() -> Option<Self> {
        which::which("pangrowth").ok().map(|path| Self { path })
    }

    /// Create from a known path.
    pub fn from_path(path: PathBuf) -> Self {
        Self { path }
    }

    /// Compute pangenome growth from a presence/absence matrix file.
    pub fn compute_growth(&self, pa_matrix: &Path) -> Result<PangrowthResult> {
        let output = Command::new(&self.path)
            .arg("growth")
            .arg("-p")
            .arg(pa_matrix)
            .output()
            .map_err(|e| Error::ExternalTool(format!("pangrowth: {}", e)))?;

        if !output.status.success() {
            return Err(Error::ExternalTool(format!(
                "pangrowth failed: {}",
                String::from_utf8_lossy(&output.stderr)
            )));
        }

        // Parse stdout for alpha, kappa, growth curve
        self.parse_growth_output(&String::from_utf8_lossy(&output.stdout))
    }

    /// Compute core genome size curve.
    pub fn compute_core(&self, pa_matrix: &Path) -> Result<Vec<(usize, f64)>> {
        let output = Command::new(&self.path)
            .arg("core")
            .arg("-p")
            .arg(pa_matrix)
            .output()
            .map_err(|e| Error::ExternalTool(format!("pangrowth: {}", e)))?;

        if !output.status.success() {
            return Err(Error::ExternalTool(format!(
                "pangrowth core failed: {}",
                String::from_utf8_lossy(&output.stderr)
            )));
        }

        self.parse_core_output(&String::from_utf8_lossy(&output.stdout))
    }

    fn parse_growth_output(&self, output: &str) -> Result<PangrowthResult> {
        // Parse Heaps' law fit: alpha and kappa from output
        // Parse growth curve data: n pangenome_size
        let mut alpha = 0.0f64;
        let mut kappa = 0.0f64;
        let mut growth_curve = Vec::new();

        for line in output.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            // Parse based on pangrowth output format
            // Growth output has columns: n expected_pangenome_size
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                if let (Ok(n), Ok(size)) = (parts[0].parse::<usize>(), parts[1].parse::<f64>()) {
                    growth_curve.push((n, size));
                }
            }
        }

        let classification = if alpha > 0.0 {
            OpennessClassification::Open
        } else {
            OpennessClassification::Closed
        };

        Ok(PangrowthResult {
            alpha, kappa, classification, growth_curve, core_curve: vec![],
        })
    }

    fn parse_core_output(&self, output: &str) -> Result<Vec<(usize, f64)>> {
        let mut core_curve = Vec::new();
        for line in output.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                if let (Ok(n), Ok(size)) = (parts[0].parse::<usize>(), parts[1].parse::<f64>()) {
                    core_curve.push((n, size));
                }
            }
        }
        Ok(core_curve)
    }
}
```

- [ ] **Step 2: Register module and add re-exports**

In `src/downstream/evolution/mod.rs`:

```rust
pub mod pangrowth;
```

In `src/downstream/mod.rs`, add:

```rust
pub use evolution::pangrowth::{PangrowthRunner, PangrowthResult, OpennessClassification};
```

- [ ] **Step 3: Add `--pangrowth` flag to analyze subcommand**

In `src/main.rs`, add to the `Analyze` variant:

```rust
/// Run Pangrowth pangenome openness estimation
#[arg(long)]
pangrowth: bool,
```

Wire in the handler:

```rust
if analyze_args.pangrowth {
    let runner = PangrowthRunner::detect()
        .ok_or_else(|| Error::ExternalTool("pangrowth not found".into()))?;
    let pa_matrix = analyze_args.input.join("gene_presence_absence.Rtab");
    let result = runner.compute_growth(&pa_matrix)?;
    println!("Pangenome openness: {:?}", result.classification);
    println!("Heaps' law alpha: {:.4}", result.alpha);
    println!("Heaps' law kappa: {:.4}", result.kappa);
    // Write growth curve CSV
}
```

- [ ] **Step 4: Write test and commit**

```bash
cargo test --features full
git add src/downstream/evolution/pangrowth.rs src/downstream/evolution/mod.rs src/downstream/mod.rs src/main.rs
git commit -m "feat: add Pangrowth pangenome openness estimation subprocess

Wires pangrowth growth/core commands for exact pangenome openness estimation.
Classifies pangenomes as open (alpha > 0) or closed (alpha <= 0) using
Heaps' law fitting. Accessible via panminer analyze --pangrowth."
```

---

### Task 2: Add abundance visualization (HTML/D3.js)

**Files:**
- Create: `src/output/abundance_viz.rs`
- Modify: `src/output/mod.rs` (register module)
- Modify: `src/main.rs` (add `--abundance` flag to analyze subcommand)
- Modify: `src/config.rs` (add Abundance to OutputFormat enum)

- [ ] **Step 1: Create AbundanceVizWriter module**

Create `src/output/abundance_viz.rs`:

```rust
//! Gene frequency and rarefaction curve visualization.
//!
//! Generates an HTML report with D3.js plots:
//! - U-shape plot: gene families vs. number of genomes
//! - Rarefaction curve: cumulative pangenome size vs. genomes added
//! - Heaps' law fit overlay
//! - Core/soft-core/shell/cloud partition bars

use std::path::Path;
use crate::error::Result;
use crate::graph::BitPackedMatrix;

pub struct AbundanceVizWriter;

impl AbundanceVizWriter {
    /// Generate an HTML report with gene frequency and rarefaction plots.
    pub fn write_report(
        matrix: &BitPackedMatrix,
        output_path: &Path,
        heaps_fit: Option<&crate::downstream::exploration::accumulation::HeapsLawFit>,
    ) -> Result<()> {
        // 1. Compute gene frequency histogram from BitPackedMatrix
        // 2. Generate D3.js HTML with:
        //    a. U-shape bar chart (x: #genomes, y: #gene_families)
        //    b. Rarefaction line chart
        //    c. Heaps' law overlay if available
        //    d. Partition bars (core/soft-core/shell/cloud)
        // 3. Write HTML to output_path
        todo!("implement")
    }

    /// Compute gene frequency histogram from presence/absence matrix.
    fn compute_gene_frequency(matrix: &BitPackedMatrix) -> Vec<(usize, usize)> {
        let num_genomes = matrix.num_genomes();
        let mut freq_counts: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        for cluster_idx in 0..matrix.num_clusters() {
            let count = matrix.count_present(cluster_idx);
            *freq_counts.entry(count).or_insert(0) += 1;
        }
        let mut result: Vec<(usize, usize)> = freq_counts.into_iter().collect();
        result.sort_by_key(|(k, _)| *k);
        result
    }
}
```

- [ ] **Step 2: Generate D3.js HTML template**

Follow the same pattern as `src/output/qc_viz.rs`:
- Embed D3.js from CDN (`d3js.org/d3.v7.min.js`)
- Create SVG bar chart for gene frequency (U-shape)
- Create SVG line chart for rarefaction curve
- Add Heaps' law fit line overlay
- Add partition bars with color coding

The HTML template should include embedded data as JSON arrays and render with D3.js.

- [ ] **Step 3: Add `Abundance` to OutputFormat and `--abundance` to analyze**

In `src/config.rs`:

```rust
pub enum OutputFormat {
    // ... existing variants ...
    Abundance,
}
```

Parse `"abundance"` in `parse_formats`.

In `src/main.rs`, add to `Analyze`:

```rust
/// Generate gene abundance visualization
#[arg(long)]
abundance: bool,
```

- [ ] **Step 4: Register module and wire into output pipeline**

In `src/output/mod.rs`:

```rust
pub mod abundance_viz;
```

Add to the output pipeline:

```rust
if formats.contains(&OutputFormat::Abundance) {
    abundance_viz::AbundanceVizWriter::write_report(
        &matrix, &paths.output_dir.join("abundance_report.html"), heaps_fit.as_ref()
    )?;
}
```

- [ ] **Step 5: Write test and commit**

```bash
cargo test --features full
git add src/output/abundance_viz.rs src/output/mod.rs src/config.rs src/main.rs
git commit -m "feat: add gene abundance HTML visualization

Generates interactive D3.js reports with U-shape gene frequency plot,
rarefaction curves, Heaps' law fit overlay, and partition bars. Accessible
via --formats abundance or panminer analyze --abundance."
```

---

### Task 3: Add Orphos/Prodigal gene calling (feature-gated)

**Files:**
- Create: `src/io/orphos.rs` (feature-gated behind `prodigal`)
- Modify: `src/io/mod.rs` (register module, feature-gated)
- Modify: `Cargo.toml` (add `orphos-core` optional dependency and `prodigal` feature)
- Modify: `src/main.rs` (add `--pipeline-mode prodigal` option)
- Modify: `src/pipeline.rs` (wire OrphosRunner)

- [ ] **Step 1: Add orphos-core dependency to Cargo.toml**

In `Cargo.toml`, add:

```toml
[dependencies]
orphos-core = { version = "0.1.0", optional = true }

[features]
prodigal = ["orphos-core"]
```

- [ ] **Step 2: Create OrphosRunner module**

Create `src/io/orphos.rs`:

```rust
//! Orphos/Prodigal gene calling for unannotated genome assemblies.
//!
//! Uses the orphos-core Rust crate (feature-gated) to predict
//! protein-coding genes in bacterial/archaeal genomes without
//! requiring external Prokka annotation.

use crate::error::Result;

/// Predicted gene from Orphos.
pub struct PredictedGene {
    pub gene_id: String,
    pub contig: String,
    pub start: usize,
    pub end: usize,
    pub strand: crate::graph::Strand,
    pub sequence: Vec<u8>,
    pub protein: Vec<u8>,
}

pub struct OrphosRunner {
    metagenomic: bool,
    closed_ends: bool,
}

impl OrphosRunner {
    pub fn new() -> Self {
        Self {
            metagenomic: false,
            closed_ends: false,
        }
    }

    pub fn with_metagenomic(mut self, meta: bool) -> Self {
        self.metagenomic = meta;
        self
    }

    pub fn with_closed_ends(mut self, closed: bool) -> Self {
        self.closed_ends = closed;
        self
    }

    /// Predict genes in a FASTA file.
    pub fn predict_genes(&self, fasta_path: &std::path::Path) -> Result<Vec<PredictedGene>> {
        #[cfg(feature = "prodigal")]
        {
            use orphos_core::OrphosAnalyzer;
            use orphos_core::config::OrphosConfig;

            let config = OrphosConfig {
                metagenomic: self.metagenomic,
                closed_ends: self.closed_ends,
                ..Default::default()
            };
            let mut analyzer = OrphosAnalyzer::new(config);
            let results = analyzer.analyze_file(fasta_path)
                .map_err(|e| crate::error::Error::ExternalTool(format!("orphos: {}", e)))?;

            Ok(results.genes.into_iter().map(|g| PredictedGene {
                gene_id: g.id,
                contig: g.contig.clone(),
                start: g.start,
                end: g.end,
                strand: if g.strand == orphos_core::types::Strand::Forward {
                    crate::graph::Strand::Plus
                } else {
                    crate::graph::Strand::Minus
                },
                sequence: g.dna_sequence,
                protein: g.protein_sequence,
            }).collect())
        }

        #[cfg(not(feature = "prodigal"))]
        {
            Err(crate::error::Error::FeatureNotEnabled(
                "Prodigal/Orphos gene calling requires the 'prodigal' feature flag. \
                 Rebuild with --features prodigal".into()
            ))
        }
    }

    /// Check if Orphos is available.
    pub fn is_available() -> bool {
        cfg!(feature = "prodigal")
    }
}
```

- [ ] **Step 3: Register module (feature-gated)**

In `src/io/mod.rs`:

```rust
#[cfg(feature = "prodigal")]
pub mod orphos;
```

In `src/lib.rs`:

```rust
#[cfg(feature = "prodigal")]
pub use io::orphos::OrphosRunner;
```

- [ ] **Step 4: Add `--pipeline-mode prodigal` to CLI**

In `src/main.rs`, update the `--pipeline-mode` argument to accept `prodigal`:

```rust
/// Pipeline mode: gff, dbg, or prodigal
#[arg(long, default_value = "gff")]
pipeline_mode: String,
```

In the pipeline, handle `prodigal` mode:

```rust
match pipeline_mode.as_str() {
    "gff" => PipelineMode::Gff,
    "dbg" => PipelineMode::Dbg,
    "prodigal" => PipelineMode::Prodigal,
    _ => PipelineMode::Gff,
}
```

- [ ] **Step 5: Add Prodigal pipeline mode**

In `src/config.rs`:

```rust
pub enum PipelineMode {
    Gff,
    Dbg,
    Prodigal,
}
```

In `src/pipeline.rs`, handle the `Prodigal` mode:

```rust
PipelineMode::Prodigal => {
    // 1. Read input FASTA files
    // 2. Run OrphosRunner::predict_genes on each
    // 3. Convert PredictedGene to Gene structs
    // 4. Continue with clustering and graph construction
}
```

- [ ] **Step 6: Write test (feature-gated)**

```rust
#[cfg(feature = "prodigal")]
#[test]
fn test_orphos_predicts_genes() {
    // Create a small FASTA file with a known sequence
    // Run OrphosRunner::predict_genes
    // Verify predicted genes have reasonable coordinates
}
```

- [ ] **Step 7: Run tests and commit**

```bash
cargo test --features full
git add src/io/orphos.rs src/io/mod.rs src/lib.rs Cargo.toml src/config.rs src/pipeline.rs src/main.rs
git commit -m "feat: add Orphos/Prodigal gene calling (feature-gated)

Adds orphos-core as an optional Rust dependency behind the 'prodigal' feature
flag. Enables --pipeline-mode prodigal for gene prediction on unannotated
FASTA assemblies without requiring external Prokka annotation."
```