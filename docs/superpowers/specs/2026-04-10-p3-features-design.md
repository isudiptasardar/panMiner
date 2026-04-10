# PanMiner P3 Features Design

> **Date**: 2026-04-10  
> **Author**: Implementation planning session  
> **Status**: Approved for implementation

---

## Overview

This document describes the implementation plan for PanMiner P3 (Nice-to-Have) features to achieve full Panaroo feature parity and production readiness:

1. **SIMD sequence comparison** - Replace scalar loops with actual SIMD intrinsics (AVX2/NEON)
2. **Evolutionary models (IMG/FMG)** - Gene gain/loss rate estimation
3. **Mash distance estimation** - K-mer sketching for contamination detection
4. **Scoary/SpydrPick integration** - Additional GWAS/epistasis tools
5. **Container/Docker support** - Reproducibility with proper packaging
6. **Python API/PyO3 bindings** - Full bindings (not just stub)
7. **Real MSA output** - Actually invoke MAFFT/PRANK/Clustal for alignments

---

## 1. SIMD Sequence Comparison

### 1.1 Feature Description

PanMiner currently uses scalar loops for sequence comparison in fragment merging. Replacing these with SIMD intrinsics (AVX2 for x86, NEON for ARM) will provide significant performance improvements.

### 1.2 Algorithm

Use `std::intrinsics::simd` or `std::arch` module for target-specific SIMD operations:

```rust
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

fn compare_sequences_simd(a: &[u8], b: &[u8]) -> f64 {
    // Process 32 bytes at a time using AVX2
    // Or 16 bytes at a time using NEON
    let mut matches = 0;
    let mut i = 0;

    unsafe {
        while i + 31 < a.len() && i + 31 < b.len() {
            // Load 32 bytes from each sequence
            let va = _mm256_loadu_si256(a.as_ptr().add(i) as *const _);
            let vb = _mm256_loadu_si256(b.as_ptr().add(i) as *const _);
            
            // Compare and count matches
            let eq = _mm256_cmpeq_epi8(va, vb);
            // ... popcount implementation
            matches += _mm256_movemask_epi8(eq) as f64;
            i += 32;
        }
    }

    // Handle remaining bytes with scalar
    for j in i..a.len().min(b.len()) {
        if a[j] == b[j] { matches += 1.0; }
    }

    matches / a.len().max(b.len()) as f64
}
```

### 1.3 Implementation Files

- Create: `src/correction/simd.rs` - SIMD utilities for sequence comparison
- Modify: `src/correction/fragment.rs` - Use SIMD for identity calculation
- Modify: `src/correction/missing.rs` - Use SIMD for k-mer matching

### 1.4 Testing Strategy

- Unit tests comparing SIMD vs scalar results
- Benchmark tests measuring performance improvement
- Test on different architectures (x86_64, aarch64)

---

## 2. Evolutionary Models (IMG/FMG)

### 2.1 Feature Description

Implement evolutionary models for gene gain/loss rate estimation, similar to Panaroo's IMG (Infinity Maximum Pseudolikelihood for Gene family evolution) and FMG (Fixation Model of Gene family evolution).

### 2.2 Algorithm

Implement a simplified version of Panaroo's evolutionary models:

```rust
/// Evolutionary model for gene family evolution
pub struct EvolutionaryModel {
    pub name: String,
    pub gain_rate: f64,
    pub loss_rate: f64,
    pub present_species: Vec<String>,
}

impl EvolutionaryModel {
    /// Calculate likelihood of observed presence/absence pattern
    pub fn likelihood(&self, pattern: &[bool]) -> f64 { /* ... */ }
    
    /// Estimate gain/loss rates from data
    pub fn estimate_rates(&self, matrix: &BitPackedMatrix) -> (f64, f64) { /* ... */ }
}
```

### 2.3 Implementation Files

- Create: `src/evolution/mod.rs` - Main module
- Create: `src/evolution/img.rs` - IMG model implementation
- Create: `src/evolution/fmg.rs` - FMG model implementation
- Create: `src/evolution/estimator.rs` - Rate estimation logic

### 2.4 Output Format

```
# Evolutionary Model Output
# Model: IMG/FMG
# Gain Rate: 0.05
# Loss Rate: 0.02
# LogLikelihood: -123.45
cluster_id	gain_rate	loss_rate	log_likelihood
cluster_0001	0.048	0.021	-123.45
cluster_0002	0.052	0.019	-145.67
```

### 2.5 Testing Strategy

- Unit tests for likelihood calculations
- Test against known evolutionary scenarios
- Compare with Panaroo results (when available)

---

## 3. Mash Distance Estimation

### 3.1 Feature Description

Implement Mash distance estimation using MinHash sketches for contamination detection and MDS projection, similar to Panaroo's Mash-based QC.

### 3.2 Algorithm

```rust
/// MinHash sketch for k-mer set
pub struct MinHashSketch {
    hash_values: Vec<u64>,
    k: usize,
    num_hashes: usize,
}

impl MinHashSketch {
    pub fn new(k: usize, num_hashes: usize) -> Self { /* ... */ }
    
    pub fn from_sequences(sequences: &[&[u8]], k: usize, num_hashes: usize) -> Self { /* ... */ }
    
    /// Calculate Mash distance between two sketches
    pub fn mash_distance(&self, other: &MinHashSketch) -> f64 {
        let shared = self.hash_values.iter()
            .zip(&other.hash_values)
            .filter(|(a, b)| a == b)
            .count();
        1.0 - (shared as f64 / self.hash_values.len() as f64)
    }
}
```

### 3.3 Implementation Files

- Create: `src/io/mash.rs` - MinHash sketch implementation
- Create: `src/io/mash_qc.rs` - QC integration
- Modify: `src/io/qc_traits.rs` - Add Mash QC trait
- Modify: `src/pipeline.rs` - Add Mash QC to pipeline

### 3.4 Output Files

- `mash_distances.tsv` - Pairwise Mash distances
- `contamination_plot.png` - MDS projection plot (PNG/SVG)
- `contamination_bar.png` - Contamination bar chart

### 3.5 Testing Strategy

- Unit tests for MinHash sketch construction
- Test Mash distance calculations
- Test contamination detection threshold

---

## 4. Scoary Integration

### 4.1 Feature Description

Integrate Scoary for gene-phenotype association testing, similar to Panaroo's Scoary wrapper.

### 4.2 Implementation Files

- Create: `src/gwas/scoary.rs` - Scoary wrapper
- Modify: `src/gwas/mod.rs` - Add Scoary to GWAS runner enum

### 4.3 Algorithm

```rust
pub struct ScoaryRunner {
    phenotypes: Option<PathBuf>,
    random_effect: Option<String>,
}

impl GWASRunner for ScoaryRunner {
    fn run(&self, graph: &PangenomeGraph, matrix: &BitPackedMatrix) -> Result<GWASOutput> {
        // Generate input files from graph/matrix
        // Call scoary command
        // Parse output
    }
}
```

### 4.4 Testing Strategy

- Test file generation from graph/matrix
- Test command invocation (when Scoary available)
- Parse Scoary output format

---

## 5. SpydrPick Integration

### 5.1 Feature Description

Integrate SpydrPick for epistasis detection (correlated gene presence/absence).

### 5.2 Implementation Files

- Create: `src/gwas/spydrpick.rs` - SpydrPick wrapper
- Modify: `src/gwas/mod.rs` - Add SpydrPick to GWAS runner enum

### 5.3 Algorithm

SpydrPick identifies correlated gene pairs based on presence/absence patterns across genomes.

```rust
pub struct SpydrPickRunner {
    min_correlation: f64,
    p_value_threshold: f64,
}

impl SpydrPickRunner {
    pub fn find_epistatic_pairs(&self, matrix: &BitPackedMatrix) -> Vec<(String, String, f64)> {
        // Calculate correlation matrix
        // Filter by significance
        // Return correlated pairs
    }
}
```

### 5.4 Testing Strategy

- Test correlation matrix calculation
- Test epistatic pair filtering
- Verify against known epistatic pairs

---

## 6. Container/Docker Support

### 6.1 Feature Description

Create Docker/Singularity container definitions for reproducibility.

### 6.2 Implementation Files

- Create: `Dockerfile` - Docker container definition
- Create: `singularity.def` - Singularity container definition
- Create: `container/` directory for container-related files
- Modify: `README.md` - Add container instructions

### 6.3 Dockerfile Structure

```dockerfile
FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    rustc \
    cargo \
    mmseqs2 \
    mafft \
    prank \
    clustalw \
    && rm -rf /var/lib/apt/lists/*

# Install CheckM2
RUN conda install -c bioconda checkm2

# Copy and build PanMiner
WORKDIR /app
COPY . /app/
RUN cargo build --release

# Set entrypoint
ENTRYPOINT ["/app/target/release/panminer"]
```

### 6.4 Testing Strategy

- Build Docker image and verify it runs
- Test pipeline execution in container
- Verify all output formats work

---

## 7. Python API/PyO3 Bindings

### 7.1 Feature Description

Implement full PyO3 bindings for PanMiner, making it available as a Python library.

### 7.2 Implementation Files

- Create: `src/python/mod.rs` - Python module entry point
- Create: `src/python/wrapper.rs` - Type wrappers
- Modify: `Cargo.toml` - Add pyo3 dependencies
- Create: `python/panminer/` - Python package structure

### 7.3 Python API Structure

```python
import panminer

# Configure pipeline
config = panminer.PanminerConfig() \
    .with_input_files(["genome1.gff", "genome2.gff"]) \
    .with_output_dir("output") \
    .with_threads(8)

# Run pipeline
pipeline = panminer.PanminerPipeline(config)
result = pipeline.run()

# Access results
print(f"Output: {result.output_dir}")
print(f"Matrix shape: {result.matrix.shape}")
```

### 7.4 Testing Strategy

- Test Python package installation
- Test basic pipeline execution from Python
- Test type conversions (Rust ↔ Python)
- Test error handling across FFI boundary

---

## 8. Real MSA Output

### 8.1 Feature Description

Actually invoke MAFFT/PRANK/Clustal for real multiple sequence alignment output instead of placeholder metadata.

### 8.2 Implementation Files

- Modify: `src/output/alignment.rs` - Add subprocess invocation
- Modify: `src/clustering/alignment_traits.rs` - Add proper trait methods
- Modify: `src/pipeline.rs` - Add alignment step to pipeline

### 8.3 Algorithm

```rust
pub struct AlignmentWriter {
    tool: AlignmentTool,
}

impl AlignmentWriter {
    pub fn write_core_alignment(&self, graph: &PangenomeGraph, path: &Path) -> Result<()> {
        // Extract core gene sequences
        // Invoke MAFFT/PRANK/Clustal via subprocess
        // Write alignment output
    }
}
```

### 8.4 Testing Strategy

- Test subprocess invocation for each alignment tool
- Verify alignment output format (FASTA, PHYLIP, etc.)
- Test with different alignment parameters

---

## Implementation Plan

### Phase 1: Core Optimizations (SIMD, Real MSA)
1. SIMD sequence comparison module
2. Real MSA subprocess invocation
3. Tests and benchmarks

### Phase 2: Downstream Analysis (GWAS Extensions)
4. Scoary integration
5. SpydrPick integration
6. GWAS result parsing

### Phase 3: Evolutionary Analysis
7. Evolutionary models module
8. Rate estimation implementation
9. Output generation

### Phase 4: Pre-processing QC
10. Mash distance estimation
11. Contamination detection
12. MDS projection

### Phase 5: Python Integration
13. PyO3 bindings
14. Python package structure
15. Python tests

### Phase 6: Reproducibility
16. Docker/Singularity definitions
17. Container testing

---

## Success Criteria

- All P3 features implemented and passing tests
- SIMD implementation shows 2x+ speedup over scalar
- MSA output matches MAFFT/PRANK/Clustal output
- Python API works for basic pipeline execution
- Docker container builds and runs successfully
- Documentation updated for all new features

---

## Notes

- SIMD implementations should gracefully degrade to scalar on unsupported architectures
- Evolutionary models may need simplification from Panaroo's full implementation
- Container should use multi-stage builds for smaller image sizes
- Python API should handle all error cases with proper Python exceptions
