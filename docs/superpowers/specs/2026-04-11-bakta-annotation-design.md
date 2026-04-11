# Bakta Annotation Integration Design

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Bakta as a pre-annotation wrapper in PanMiner, mirroring Panaroo's Prokka integration, so users can provide raw genome assemblies and get consistent annotations before pangenome analysis.

**Architecture:** Subprocess runner pattern (matching existing CheckmQcRunner/MMseqsRunner/MafftRunner). Bakta runs as an external Python tool, producing GFF3 output that PanMiner's existing GFF parser consumes. Pipeline adds an optional Phase 0.5 between QC and parsing.

**Tech Stack:** Rust (subprocess via `std::process::Command`), Bakta CLI (Python, installed externally), existing GffParser for output consumption.

---

## Context

Panaroo re-annotates all input genomes with Prokka to ensure consistent annotations before pangenome analysis. PanMiner currently only accepts pre-annotated GFF3 files, requiring users to run annotation tools themselves. Adding Bakta as a pre-annotation wrapper closes this usability gap with a modern, more accurate annotator.

Bakta is a Python-based bacterial genome annotator that produces GFF3 output compatible with PanMiner's existing parser. It uses MD5 hash-based protein identification (avoiding expensive homology searches for known proteins), provides dbxref-rich annotations, and supports both full and light database modes.

## Key Design Decisions

1. **Subprocess runner, not Python API** — Matches existing patterns (CheckM2, MMseqs2, MAFFT). No PyO3 dependency. Bakta stays a Python tool managed externally.

2. **Graceful fallback** — If Bakta is not installed, PanMiner falls back to direct GFF input. No hard requirement.

3. **Auto-download database** — Resolves DB path from flag > env var > default location, auto-downloads if missing (unless `--no-bakta-db-download`).

4. **Phase 0.5 pipeline placement** — Between QC (Phase 0) and GFF parsing (Phase 1). Re-annotation happens before any pangenome processing.

5. **Mixed input support** — GFF files pass through unchanged; FASTA/GenBank files get annotated by Bakta.

---

## Module Structure

### New Files

- `src/io/bakta.rs` — `BaktaRunner` struct: detect, resolve_db, download_db, annotate, GenBank-to-FASTA converter
- Update `src/io/mod.rs` — export `BaktaRunner`
- Update `src/config.rs` — add re-annotation config fields
- Update `src/main.rs` — add CLI flags
- Update `src/pipeline.rs` — add Phase 0.5

### BaktaRunner API

```rust
pub struct BaktaRunner {
    bakta_path: PathBuf,      // Path to bakta binary
    db_path: PathBuf,         // Path to Bakta database
    threads: usize,           // Thread count
    output_dir: PathBuf,     // Temporary output directory
    keep_contig_headers: bool, // --keep-contig-headers flag
}

impl BaktaRunner {
    pub fn detect() -> Option<Self>
    pub fn resolve_db(no_download: bool, db_type: BaktaDbType) -> Result<PathBuf>
    pub fn download_db(output: &Path, db_type: BaktaDbType) -> Result<PathBuf>
    pub fn annotate(&self, input: &Path) -> Result<PathBuf>
    pub fn annotate_batch(&self, inputs: &[PathBuf]) -> Result<Vec<PathBuf>>
    pub fn name(&self) -> &str
}

pub enum BaktaDbType { Full, Light }

fn genbank_to_fasta(input: &Path) -> Result<PathBuf>
```

### CLI Flags

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--reannotate` | `-r` | false | Re-annotate input genomes with Bakta before analysis |
| `--bakta-db` | | auto | Path to Bakta database directory |
| `--bakta-db-type` | | `full` | Database type for auto-download: `full` or `light` |
| `--bakta-threads` | | auto | Number of threads for Bakta (default: same as pipeline) |
| `--no-bakta-db-download` | | false | Fail if Bakta DB not found instead of auto-downloading |
| `--keep-bakta-output` | | false | Keep Bakta output files after pipeline completes |

### Config Changes

```rust
// In PanminerConfig:
pub reannotate: bool,              // --reannotate
pub bakta_db_path: Option<PathBuf>, // --bakta-db
pub bakta_db_type: BaktaDbType,    // --bakta-db-type
pub bakta_threads: usize,          // --bakta-threads
pub no_bakta_db_download: bool,    // --no-bakta-db-download
pub keep_bakta_output: bool,       // --keep-bakta-output
```

### Pipeline Flow

```
Phase 0:   QC (CheckM2, existing)
Phase 0.5: Re-annotation (NEW, optional, only if --reannotate)
            ├── Detect Bakta installation
            ├── Resolve/download database
            ├── For each input file:
            │   ├── .gff/.gff3 → pass through unchanged
            │   ├── .fasta/.fna/.fa → run Bakta → get .gff3
            │   └── .gbk/.gb/.gbff → convert to FASTA → run Bakta → get .gff3
            └── Collect all GFF3 paths
Phase 1:    Parse GFF3 files (existing, uses output from Phase 0.5)
Phase 2-6:  (unchanged)
```

### Error Handling

| Scenario | Behavior |
|----------|----------|
| `--reannotate` set, Bakta not installed | Log warning, skip re-annotation, use input files directly |
| `--reannotate` set, Bakta found, DB not found | Auto-download DB unless `--no-bakta-db-download` set |
| `--no-bakta-db-download` set, DB not found | Error with installation instructions |
| Bakta annotation fails for one genome | Log error, skip that genome, continue with remaining |
| Bakta annotation fails for ALL genomes | Return error |
| Mixed input (.gff + .fasta) | Annotate only FASTA files, pass GFF files through |
| `.gbk` input without Bakta | Error — GenBank files require Bakta for conversion |

### GenBank-to-FASTA Conversion

Simple parser (~50 lines) that extracts sequence data from GenBank format:
- Find `ORIGIN` line
- Read lines until `//`
- Strip line numbers, spaces, and newlines
- Write to temporary FASTA file

### Temporary File Management

- Bakta output goes to `<output_dir>/bakta_tmp/`
- Cleaned up after GFF parsing unless `--keep-bakta-output` is set
- GenBank-to-FASTA conversions go to `<output_dir>/bakta_tmp/converted/`

---

## Testing Strategy

### Unit Tests (feature-gated)

- `test_bakta_detect_found` / `test_bakta_detect_not_found` — mock `bakta --version`
- `test_bakta_resolve_db_from_flag` / `test_bakta_resolve_db_from_env` — priority order
- `test_bakta_annotate_mock` — mock subprocess producing known GFF3
- `test_genbank_to_fasta` — small `.gbk` snippet conversion

### Integration Tests

- `test_bakta_reannotation_with_mock` — temp FASTA → mock Bakta → parse GFF3
- `test_bakta_graceful_fallback` — pipeline works without Bakta
- `test_mixed_input_gff_and_fasta` — GFF pass-through + FASTA annotation
- `test_bakta_db_type_flag` — verify light/full DB type propagates

### Test Infrastructure

- Mock Bakta via shell script producing known GFF3 output (no Python dependency)
- `#[cfg(feature = "bakta")]` for tests requiring the binary
- `tempfile::tempdir()` for all output, matching existing test patterns
- Separate `tests/bakta_integration_test.rs` for optional live-Bakta CI tests

---

## Out of Scope (YAGNI)

- PyO3 integration (subprocess pattern is sufficient and matches existing tools)
- Bakta JSON output parsing (GFF3 is sufficient, JSON is internal format)
- Bakta-specific attribute extraction (`db_xref`, `locus_tag`, `ec_number`) — current parser handles basic annotation; enhanced attributes can be added later
- ProteinFASTA-only input (Bakta requires nucleotide FASTA)
- `bakta_proteins` command (not needed for pangenome analysis)
- Parallel Bakta runs across genomes (Bakta already uses `--threads`; adding Rayon parallelism on top would over-subscribe CPUs)