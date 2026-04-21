//! Bakta annotation runner for genome re-annotation.
//!
//! Runs Bakta as a subprocess to annotate raw genome assemblies (FASTA/GenBank)
//! before pangenome analysis. This mirrors Panaroo's Prokka re-annotation step
//! but uses Bakta, which provides more accurate annotations via MD5-based protein
//! identification and dbxref-rich output.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::Result;

/// Bakta database type for auto-download.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaktaDbType {
    /// Full database (~6GB, best results)
    Full,
    /// Light database (~500MB, faster runtime, fewer features)
    Light,
}

impl std::fmt::Display for BaktaDbType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BaktaDbType::Full => write!(f, "full"),
            BaktaDbType::Light => write!(f, "light"),
        }
    }
}

/// Bakta annotation runner.
///
/// Detects Bakta installation, resolves/downloads the database,
/// and runs `bakta` CLI to produce GFF3 output from raw genome assemblies.
pub struct BaktaRunner {
    /// Path to bakta binary
    bakta_path: PathBuf,
    /// Path to Bakta database directory
    db_path: PathBuf,
    /// Number of threads for Bakta (0 = auto-detect)
    threads: usize,
    /// Output directory for Bakta results
    output_dir: PathBuf,
    /// Whether to keep original contig headers
    keep_contig_headers: bool,
}

impl BaktaRunner {
    /// Create a new BaktaRunner with explicit paths.
    pub fn new(bakta_path: PathBuf, db_path: PathBuf) -> Self {
        Self {
            bakta_path,
            db_path,
            threads: 0,
            output_dir: PathBuf::from("."),
            keep_contig_headers: true,
        }
    }

    /// Set the number of threads (0 = auto-detect).
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set the output directory for Bakta results.
    pub fn with_output_dir(mut self, dir: PathBuf) -> Self {
        self.output_dir = dir;
        self
    }

    /// Set whether to keep original contig headers.
    pub fn with_keep_contig_headers(mut self, keep: bool) -> Self {
        self.keep_contig_headers = keep;
        self
    }

    /// Detect if Bakta is installed on the system.
    ///
    /// Returns `Some(BaktaRunner)` if `bakta --version` succeeds, `None` otherwise.
    pub fn detect() -> Option<Self> {
        let bakta_path = which_bakta()?;

        // Verify it actually works
        let output = Command::new(&bakta_path)
            .arg("--version")
            .output()
            .ok()?;

        if !output.status.success() {
            return None;
        }

        // Find database path
        let db_path = resolve_db_path(None);

        Some(Self {
            bakta_path,
            db_path,
            threads: 0,
            output_dir: PathBuf::from("."),
            keep_contig_headers: true,
        })
    }

    /// Resolve the Bakta database path.
    ///
    /// Priority: explicit_path > BAKTA_DB env > ~/.bakta/db
    pub fn resolve_db(explicit_path: Option<&Path>) -> PathBuf {
        resolve_db_path(explicit_path)
    }

    /// Download the Bakta database if not already present.
    ///
    /// Uses `bakta_db download` to fetch the database.
    pub fn download_db(output_path: &Path, db_type: BaktaDbType) -> Result<PathBuf> {
        let bakta_db_cmd = which_bakta_db().ok_or_else(|| {
            crate::Error::ExternalTool(
                "bakta_db command not found. Install Bakta first: `conda install -c conda-forge -c bioconda bakta`".to_string()
            )
        })?;

        let db_dir = output_path.join(format!("db-{}", db_type));

        if db_dir.exists() {
            tracing::info!("Bakta database already exists at {:?}", db_dir);
            return Ok(db_dir);
        }

        tracing::info!("Downloading Bakta {} database to {:?}", db_type, output_path);

        let output = Command::new(&bakta_db_cmd)
            .arg("download")
            .arg("--output")
            .arg(output_path)
            .arg("--type")
            .arg(db_type.to_string())
            .output()
            .map_err(|e| crate::Error::ExternalTool(format!(
                "Failed to run bakta_db download: {}", e
            )))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::ExternalTool(format!(
                "bakta_db download failed: {}", stderr.trim()
            )));
        }

        tracing::info!("Bakta database downloaded to {:?}", db_dir);
        Ok(db_dir)
    }

    /// Annotate a single genome assembly with Bakta.
    ///
    /// Takes a FASTA/GenBank file path, runs Bakta, and returns
    /// the path to the resulting GFF3 file.
    pub fn annotate(&self, input_path: &Path) -> Result<PathBuf> {
        let genome_name = input_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("genome");

        let output_dir = self.output_dir.join("bakta_tmp").join(genome_name);
        std::fs::create_dir_all(&output_dir)?;

        tracing::info!("Annotating {:?} with Bakta", input_path);

        let mut cmd = Command::new(&self.bakta_path);
        cmd.arg("--db")
            .arg(&self.db_path)
            .arg("--output")
            .arg(&output_dir)
            .arg("--skip-plot"); // Skip circular genome plots

        if self.threads > 0 {
            cmd.arg("--threads").arg(self.threads.to_string());
        }

        if self.keep_contig_headers {
            cmd.arg("--keep-contig-headers");
        }

        cmd.arg(input_path);

        let output = cmd.output().map_err(|e| {
            crate::Error::ExternalTool(format!("Failed to run Bakta: {}", e))
        })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::bakta_annotation_failed(
                genome_name, &stderr
            ));
        }

        // Find the GFF3 output file
        let gff_path = output_dir.join(format!("{}.gff3", genome_name));
        if gff_path.exists() {
            tracing::info!("Bakta annotation complete: {:?}", gff_path);
            return Ok(gff_path);
        }

        // Try to find any .gff3 file in the output directory
        if let Ok(entries) = std::fs::read_dir(&output_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if let Some(ext) = path.extension() {
                    if ext == "gff3" {
                        tracing::info!("Bakta annotation complete: {:?}", path);
                        return Ok(path);
                    }
                }
            }
        }

        Err(crate::Error::ExternalTool(format!(
            "Bakta output GFF3 not found in {:?}", output_dir
        )))
    }

    /// Annotate multiple genome assemblies in sequence.
    ///
    /// Processes each input file, collecting GFF3 output paths.
    /// Files that are already GFF3 are passed through unchanged.
    pub fn annotate_batch(&self, inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
        let mut gff_paths = Vec::new();

        for input in inputs {
            if is_gff_file(input) {
                // Pass GFF files through unchanged
                tracing::info!("Passing through GFF file: {:?}", input);
                gff_paths.push(input.clone());
            } else if is_genbank_file(input) {
                // Convert GenBank to FASTA first, then annotate
                let fasta_path = genbank_to_fasta(input)?;
                let gff_path = self.annotate(&fasta_path)?;
                gff_paths.push(gff_path);
            } else {
                // Assume FASTA — annotate directly
                let gff_path = self.annotate(input)?;
                gff_paths.push(gff_path);
            }
        }

        Ok(gff_paths)
    }

    /// Get the name of this runner.
    pub fn name(&self) -> &str {
        "Bakta"
    }

    /// Get the path to the Bakta binary.
    pub fn name_path(&self) -> PathBuf {
        self.bakta_path.clone()
    }
}

/// Find the bakta executable on the system PATH.
fn which_bakta() -> Option<PathBuf> {
    // Try direct command first
    if Command::new("bakta")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
    {
        return Some(PathBuf::from("bakta"));
    }

    // Try common conda paths
    let home = std::env::var("HOME").ok()?;
    let conda_paths = [
        format!("{}/.conda/envs/base/bin/bakta", home),
        format!("{}/miniconda3/bin/bakta", home),
        format!("{}/anaconda3/bin/bakta", home),
    ];

    for path in &conda_paths {
        let pb = PathBuf::from(path);
        if pb.exists() {
            return Some(pb);
        }
    }

    None
}

/// Find the bakta_db executable on the system PATH.
fn which_bakta_db() -> Option<PathBuf> {
    if Command::new("bakta_db")
        .arg("--help")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
    {
        return Some(PathBuf::from("bakta_db"));
    }
    None
}

/// Resolve the Bakta database path.
///
/// Priority: explicit_path > BAKTA_DB env > ~/.bakta/db
fn resolve_db_path(explicit_path: Option<&Path>) -> PathBuf {
    // 1. Explicit path from CLI flag
    if let Some(path) = explicit_path {
        if path.exists() {
            return path.to_path_buf();
        }
    }

    // 2. BAKTA_DB environment variable
    if let Ok(env_path) = std::env::var("BAKTA_DB") {
        let path = PathBuf::from(&env_path);
        if path.exists() {
            return path;
        }
    }

    // 3. Default location (Unix)
    if let Ok(home) = std::env::var("HOME") {
        let default_path = PathBuf::from(home).join(".bakta").join("db");
        if default_path.exists() {
            return default_path;
        }
    }

    // 4. Windows: check USERPROFILE
    if let Ok(userprofile) = std::env::var("USERPROFILE") {
        let default_path = PathBuf::from(userprofile).join(".bakta").join("db");
        if default_path.exists() {
            return default_path;
        }
    }

    // Return best-effort default using detected HOME/USERPROFILE
    if let Ok(home) = std::env::var("HOME") {
        PathBuf::from(home).join(".bakta").join("db")
    } else if let Ok(userprofile) = std::env::var("USERPROFILE") {
        PathBuf::from(userprofile).join(".bakta").join("db")
    } else {
        PathBuf::from(".bakta_db")
    }
}

/// Check if a file has a GFF/GFF3 extension.
pub fn is_gff_file(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .map(|e| matches!(e.to_lowercase().as_str(), "gff" | "gff3"))
        .unwrap_or(false)
}

/// Check if a file has a GenBank extension.
pub fn is_genbank_file(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .map(|e| matches!(e.to_lowercase().as_str(), "gb" | "gbk" | "gbff" | "genbank"))
        .unwrap_or(false)
}

/// Convert a GenBank file to FASTA format.
///
/// Extracts the ORIGIN section from a GenBank file and converts
/// it to a simple FASTA file suitable for Bakta input.
pub fn genbank_to_fasta(input: &Path) -> Result<PathBuf> {
    let content = std::fs::read_to_string(input)
        .map_err(crate::Error::Io)?;

    let genome_name = input
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("genome");

    // Find the ORIGIN section
    let mut sequence = String::new();
    let mut in_origin = false;

    for line in content.lines() {
        if line.starts_with("ORIGIN") {
            in_origin = true;
            continue;
        }
        if in_origin {
            if line.starts_with("//") {
                break;
            }
            // GenBank ORIGIN lines have format: "        1 atcgatcg atcgatcg ..."
            // Strip the leading number (sequence position) and spaces
            let trimmed = line.trim_start();
            // Skip the sequence position number at the start
            let after_number = trimmed.find(char::is_alphabetic).unwrap_or(0);
            let seq_part = &trimmed[after_number..];
            for ch in seq_part.chars() {
                if ch.is_ascii_alphabetic() {
                    sequence.push(ch);
                }
                // Skip spaces and digits within the sequence area
            }
        }
    }

    if sequence.is_empty() {
        return Err(crate::Error::ExternalTool(format!(
            "No ORIGIN section found in GenBank file: {:?}", input
        )));
    }

    // Write FASTA file
    let output_path = input.with_extension("fna");
    let fasta_content = format!(">{}\n{}", genome_name, sequence);
    std::fs::write(&output_path, &fasta_content)
        .map_err(crate::Error::Io)?;

    tracing::info!("Converted GenBank {:?} to FASTA {:?}", input, output_path);
    Ok(output_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bakta_db_type_display() {
        assert_eq!(BaktaDbType::Full.to_string(), "full");
        assert_eq!(BaktaDbType::Light.to_string(), "light");
    }

    #[test]
    fn test_is_gff_file() {
        assert!(is_gff_file(Path::new("genome.gff")));
        assert!(is_gff_file(Path::new("genome.gff3")));
        assert!(is_gff_file(Path::new("genome.GFF3")));
        assert!(!is_gff_file(Path::new("genome.fasta")));
        assert!(!is_gff_file(Path::new("genome.gbk")));
    }

    #[test]
    fn test_is_genbank_file() {
        assert!(is_genbank_file(Path::new("genome.gb")));
        assert!(is_genbank_file(Path::new("genome.gbk")));
        assert!(is_genbank_file(Path::new("genome.gbff")));
        assert!(is_genbank_file(Path::new("genome.genbank")));
        assert!(!is_genbank_file(Path::new("genome.fasta")));
        assert!(!is_genbank_file(Path::new("genome.gff")));
    }

    #[test]
    fn test_genbank_to_fasta() {
        let temp_dir = tempfile::tempdir().unwrap();
        let gbk_path = temp_dir.path().join("test.gb");

        let gbk_content = r#"LOCUS       test                100 bp    DNA     linear   BCT 01-JAN-2024
DEFINITION  Test genome.
ACCESSION   test
VERSION     test.1
ORIGIN
        1 atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg
       61 atcgatcgat cgatcgatcg atcgatcgat cgatcgatcg
//
"#;
        std::fs::write(&gbk_path, gbk_content).unwrap();

        let fasta_path = genbank_to_fasta(&gbk_path).unwrap();
        let fasta_content = std::fs::read_to_string(&fasta_path).unwrap();

        assert!(fasta_content.starts_with(">test\n"));
        // Verify nucleotide content is extracted (no line numbers, no spaces)
        let seq = fasta_content.strip_prefix(">test\n").unwrap();
        assert!(!seq.is_empty());
        assert!(seq.chars().all(|c| "atcgATCG".contains(c)));
    }

    #[test]
    fn test_genbank_to_fasta_missing_origin() {
        let temp_dir = tempfile::tempdir().unwrap();
        let gbk_path = temp_dir.path().join("empty.gb");
        std::fs::write(&gbk_path, "LOCUS test\nDEFINITION empty\n//").unwrap();

        let result = genbank_to_fasta(&gbk_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_bakta_runner_creation() {
        let runner = BaktaRunner::new(
            PathBuf::from("bakta"),
            PathBuf::from("/path/to/db"),
        );
        assert_eq!(runner.name(), "Bakta");
        assert_eq!(runner.threads, 0);
        assert!(runner.keep_contig_headers);
    }

    #[test]
    fn test_bakta_runner_builder_pattern() {
        let runner = BaktaRunner::new(
            PathBuf::from("bakta"),
            PathBuf::from("/path/to/db"),
        )
        .with_threads(8)
        .with_keep_contig_headers(false);

        assert_eq!(runner.threads, 8);
        assert!(!runner.keep_contig_headers);
    }

    #[test]
    fn test_resolve_db_path_with_explicit() {
        let tmp = tempfile::tempdir().unwrap();
        let db_path = tmp.path().join("db");
        std::fs::create_dir_all(&db_path).unwrap();

        let result = resolve_db_path(Some(&db_path));
        assert_eq!(result, db_path);
    }
}