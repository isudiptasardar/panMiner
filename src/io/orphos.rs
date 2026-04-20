//! Prodigal gene calling for unannotated genome assemblies.
//!
//! Uses the Prodigal binary (feature-gated) to predict protein-coding
//! genes in bacterial/archaeal genomes. This allows PanMiner to process
//! raw FASTA assemblies without requiring pre-annotated GFF3 files.
//!
//! Prodigal (PROkaryotic Dynamic Programming Gene-finding ALgorithm) is
//! a fast, reliable protein-coding gene predictor for bacterial and
//! archaeal genomes. It runs as a subprocess, producing GFF, protein
//! FASTA, and nucleotide FASTA output files.
//!
//! # Feature flag
//!
//! This module requires the `prodigal` feature:
//!
//! ```toml
//! panminer = { features = ["prodigal"] }
//! ```
//!
//! # Example
//!
//! ```no_run
//! #[cfg(feature = "prodigal")]
//! {
//!     use panminer::io::OrphosRunner;
//!
//!     if OrphosRunner::is_installed() {
//!         let runner = OrphosRunner::detect().unwrap();
//!         let genes = runner.predict_genes(std::path::Path::new("genome.fna")).unwrap();
//!         println!("Predicted {} genes", genes.len());
//!     }
//! }
//! ```

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::Result;

/// Strand direction for a predicted gene.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// Forward strand (+)
    Forward,
    /// Reverse strand (-)
    Reverse,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

/// A predicted gene from Prodigal output.
#[derive(Debug, Clone)]
pub struct PredictedGene {
    /// Gene identifier (e.g., "PRODIGAL_00001" from the GFF attributes)
    pub gene_id: String,
    /// Contig / sequence name where this gene is located
    pub contig: String,
    /// 1-based start position
    pub start: usize,
    /// 1-based end position (inclusive)
    pub end: usize,
    /// Strand orientation
    pub strand: Strand,
    /// Nucleotide sequence (from the .fna output)
    pub sequence: Vec<u8>,
    /// Protein sequence (from the .faa output)
    pub protein: Vec<u8>,
}

/// Prodigal gene-calling runner.
///
/// Detects the Prodigal binary, constructs the appropriate command line,
/// runs the subprocess, and parses the output files into `PredictedGene`
/// structs.
///
/// Uses a builder pattern for configuration:
///
/// ```no_run
/// #[cfg(feature = "prodigal")]
/// {
///     use panminer::io::OrphosRunner;
///
///     let runner = OrphosRunner::new()
///         .with_metagenomic(true)
///         .with_closed_ends(false);
/// }
/// ```
#[cfg(feature = "prodigal")]
pub struct OrphosRunner {
    /// Path to the prodigal binary
    prodigal_path: PathBuf,
    /// Use metagenomic mode (-p meta)
    metagenomic: bool,
    /// Use closed ends (-c) — assumes the input sequence ends are complete
    closed_ends: bool,
}

#[cfg(feature = "prodigal")]
impl OrphosRunner {
    /// Create a new OrphosRunner with default settings.
    ///
    /// Defaults: normal mode (not metagenomic), open ends.
    /// The prodigal binary is located via `which::which`.
    pub fn new() -> Self {
        let prodigal_path = which::which("prodigal")
            .unwrap_or_else(|_| PathBuf::from("prodigal"));

        Self {
            prodigal_path,
            metagenomic: false,
            closed_ends: false,
        }
    }

    /// Detect Prodigal on the system PATH.
    ///
    /// Returns `Some(OrphosRunner)` if `prodigal --version` succeeds,
    /// `None` otherwise.
    pub fn detect() -> Option<Self> {
        let prodigal_path = which::which("prodigal").ok()?;

        // Verify it actually works by running --version
        let output = Command::new(&prodigal_path)
            .arg("--version")
            .output()
            .ok()?;

        if !output.status.success() {
            return None;
        }

        Some(Self {
            prodigal_path,
            metagenomic: false,
            closed_ends: false,
        })
    }

    /// Set metagenomic mode.
    ///
    /// When enabled, Prodigal runs in metagenomic mode (`-p meta`),
    /// which is suitable for fragmented genome assemblies,
    /// metagenomic contigs, or viral genomes where the standard
    /// training step is not possible.
    pub fn with_metagenomic(mut self, meta: bool) -> Self {
        self.metagenomic = meta;
        self
    }

    /// Set closed ends mode.
    ///
    /// When enabled, Prodigal assumes the sequence ends are complete
    /// (`-c`). This can improve predictions at sequence boundaries
    /// for complete genomes or closed contigs.
    pub fn with_closed_ends(mut self, closed: bool) -> Self {
        self.closed_ends = closed;
        self
    }

    /// Check if Prodigal is available for this runner.
    ///
    /// Returns `true` if the stored prodigal path points to an
    /// executable binary.
    pub fn is_available(&self) -> bool {
        Command::new(&self.prodigal_path)
            .arg("--version")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    }

    /// Predict genes in a FASTA assembly using Prodigal.
    ///
    /// Runs Prodigal as a subprocess, parses the GFF, protein FASTA,
    /// and nucleotide FASTA output files, and returns a vector of
    /// `PredictedGene` structs combining all information.
    ///
    /// # Arguments
    ///
    /// * `fasta_path` - Path to the input genome assembly in FASTA format
    ///
    /// # Errors
    ///
    /// Returns `Error::ExternalTool` if Prodigal fails to run or exits
    /// with a non-zero status code. Returns `Error::Io` if temp directory
    /// creation fails.
    pub fn predict_genes(&self, fasta_path: &Path) -> Result<Vec<PredictedGene>> {
        let temp_dir = tempfile::tempdir()?;
        let genome_name = fasta_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("genome");

        let gff_path = temp_dir.path().join(format!("{}.gff", genome_name));
        let protein_path = temp_dir.path().join(format!("{}.faa", genome_name));
        let nucleotide_path = temp_dir.path().join(format!("{}.fna", genome_name));

        tracing::info!("Predicting genes in {:?} with Prodigal", fasta_path);

        // Build the Prodigal command
        let mut cmd = Command::new(&self.prodigal_path);
        cmd.arg("-i").arg(fasta_path)
            .arg("-a").arg(&protein_path)
            .arg("-d").arg(&nucleotide_path)
            .arg("-o").arg(&gff_path)
            .arg("-f").arg("gff");

        if self.metagenomic {
            cmd.arg("-p").arg("meta");
        }

        if self.closed_ends {
            cmd.arg("-c");
        }

        let output = cmd.output().map_err(|e| {
            crate::Error::ExternalTool(format!("Failed to run Prodigal: {}", e))
        })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::ExternalTool(format!(
                "Prodigal gene prediction failed: {}", stderr.trim()
            )));
        }

        // Parse the output files
        let gene_coords = parse_prodigal_gff(&gff_path)?;
        let protein_seqs = parse_fasta_sequences(&protein_path)?;
        let nucleotide_seqs = parse_fasta_sequences(&nucleotide_path)?;

        // Combine into PredictedGene structs
        let mut genes = Vec::with_capacity(gene_coords.len());

        for (gene_id, coord) in gene_coords {
            let protein = protein_seqs
                .get(&gene_id)
                .cloned()
                .unwrap_or_default();
            let sequence = nucleotide_seqs
                .get(&gene_id)
                .cloned()
                .unwrap_or_default();

            genes.push(PredictedGene {
                gene_id,
                contig: coord.contig,
                start: coord.start,
                end: coord.end,
                strand: coord.strand,
                sequence,
                protein,
            });
        }

        tracing::info!("Prodigal predicted {} genes", genes.len());

        // Temp directory is cleaned up automatically when temp_dir is dropped
        Ok(genes)
    }

    /// Get the name of this runner.
    pub fn name(&self) -> &str {
        "Prodigal"
    }

    /// Get the path to the Prodigal binary.
    pub fn binary_path(&self) -> &Path {
        &self.prodigal_path
    }

    /// Static check: is Prodigal installed on the system PATH?
    pub fn is_installed() -> bool {
        which::which("prodigal").is_ok()
    }
}

#[cfg(feature = "prodigal")]
impl Default for OrphosRunner {
    fn default() -> Self {
        Self::new()
    }
}

/// Parsed gene coordinates from a Prodigal GFF file.
struct GeneCoord {
    contig: String,
    start: usize,
    end: usize,
    strand: Strand,
}

/// Parse a Prodigal GFF output file to extract gene coordinates.
///
/// Prodigal GFF rows have the format:
/// ```text
/// contig_id  Prodigal_v2.6.3  CDS  start  end  score  strand  phase  attributes
/// ```
///
/// The gene ID is extracted from the `ID=` attribute.
#[cfg(feature = "prodigal")]
fn parse_prodigal_gff(gff_path: &Path) -> Result<HashMap<String, GeneCoord>> {
    let content = std::fs::read_to_string(gff_path)
        .map_err(crate::Error::Io)?;

    let mut coords = HashMap::new();

    for line in content.lines() {
        // Skip comments and headers
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        // Only process CDS features
        if fields[2] != "CDS" {
            continue;
        }

        let contig = fields[0].to_string();

        let start: usize = match fields[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let end: usize = match fields[4].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let strand = match fields[6] {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => continue,
        };

        // Extract gene ID from attributes (ID=...;)
        let gene_id = extract_gff_attribute(fields[8], "ID")
            .unwrap_or_else(|| format!("gene_{}_{}", contig, start));

        coords.insert(gene_id, GeneCoord {
            contig,
            start,
            end,
            strand,
        });
    }

    Ok(coords)
}

/// Extract a named attribute from a GFF attribute string.
///
/// GFF attributes are semicolon-separated key=value pairs.
/// Returns `None` if the attribute is not found.
#[cfg(feature = "prodigal")]
fn extract_gff_attribute(attributes: &str, key: &str) -> Option<String> {
    for pair in attributes.split(';') {
        let pair = pair.trim();
        if let Some(eq_pos) = pair.find('=') {
            let k = &pair[..eq_pos];
            if k == key {
                return Some(pair[eq_pos + 1..].to_string());
            }
        }
    }
    None
}

/// Parse a FASTA file into a map of sequence ID to sequence bytes.
///
/// Sequences are concatenated across line boundaries (FASTA
/// line-wrapping is handled). Sequence IDs are taken from the
/// header line after the `>` character, up to the first whitespace.
#[cfg(feature = "prodigal")]
fn parse_fasta_sequences(fasta_path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    let content = std::fs::read_to_string(fasta_path)
        .map_err(crate::Error::Io)?;

    let mut sequences = HashMap::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = Vec::new();

    for line in content.lines() {
        if line.starts_with('>') {
            // Save the previous sequence
            if let Some(id) = current_id.take() {
                sequences.insert(id, std::mem::take(&mut current_seq));
            }

            // Parse new header: >ID optional_description
            let header = line[1..].trim();
            let id = header
                .split_whitespace()
                .next()
                .unwrap_or(header)
                .to_string();
            current_id = Some(id);
        } else if current_id.is_some() {
            // Append sequence line (skip whitespace)
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                current_seq.extend_from_slice(trimmed.as_bytes());
            }
        }
    }

    // Save the last sequence
    if let Some(id) = current_id {
        sequences.insert(id, current_seq);
    }

    Ok(sequences)
}

// --- Non-feature stubs ---

#[cfg(not(feature = "prodigal"))]
pub struct OrphosRunner;

#[cfg(not(feature = "prodigal"))]
impl OrphosRunner {
    /// Returns `false` — Prodigal support is not enabled.
    pub fn is_available() -> bool {
        false
    }

    /// Returns `false` — Prodigal support is not enabled.
    pub fn is_installed() -> bool {
        false
    }
}

// --- Tests ---

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "prodigal")]
    #[test]
    fn test_orphos_detect() {
        // Just test that detect() doesn't panic
        let _ = OrphosRunner::detect();
    }

    #[cfg(feature = "prodigal")]
    #[test]
    fn test_orphos_is_installed() {
        // Just test that is_installed() doesn't panic
        let _ = OrphosRunner::is_installed();
    }

    #[cfg(feature = "prodigal")]
    #[test]
    fn test_orphos_new() {
        let runner = OrphosRunner::new();
        assert_eq!(runner.name(), "Prodigal");
        assert!(!runner.metagenomic);
        assert!(!runner.closed_ends);
    }

    #[cfg(feature = "prodigal")]
    #[test]
    fn test_orphos_builder_pattern() {
        let runner = OrphosRunner::new()
            .with_metagenomic(true)
            .with_closed_ends(true);
        assert!(runner.metagenomic);
        assert!(runner.closed_ends);
    }

    #[cfg(feature = "prodigal")]
    #[test]
    fn test_strand_display() {
        assert_eq!(format!("{}", Strand::Forward), "+");
        assert_eq!(format!("{}", Strand::Reverse), "-");
    }

    #[cfg(feature = "prodigal")]
    #[test]
    fn test_extract_gff_attribute() {
        let attrs = "ID=PRODIGAL_00001;partial=00;start_type=ATG;rbs_motif=GGAG";
        assert_eq!(
            extract_gff_attribute(attrs, "ID"),
            Some("PRODIGAL_00001".to_string())
        );
        assert_eq!(
            extract_gff_attribute(attrs, "partial"),
            Some("00".to_string())
        );
        assert_eq!(
            extract_gff_attribute(attrs, "start_type"),
            Some("ATG".to_string())
        );
        assert_eq!(extract_gff_attribute(attrs, "missing"), None);
    }

    #[cfg(feature = "prodigal")]
    #[test]
    fn test_parse_fasta_sequences() {
        let temp_dir = tempfile::tempdir().unwrap();
        let fasta_path = temp_dir.path().join("test.fna");

        let content = ">gene_1\nATGCGTACGT\nTTGCGT\n>gene_2\nGGGAAA\n";
        std::fs::write(&fasta_path, content).unwrap();

        let seqs = parse_fasta_sequences(&fasta_path).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(
            seqs.get("gene_1").unwrap(),
            &b"ATGCGTACGTTTGCGT".to_vec()
        );
        assert_eq!(
            seqs.get("gene_2").unwrap(),
            &b"GGGAAA".to_vec()
        );
    }

    #[cfg(feature = "prodigal")]
    #[test]
    fn test_parse_prodigal_gff() {
        let temp_dir = tempfile::tempdir().unwrap();
        let gff_path = temp_dir.path().join("test.gff");

        let content = r#"##gff-version 3
contig_1	Prodigal_v2.6.3	CDS	101	300	50.2	+	0	ID=PRODIGAL_00001;partial=00;start_type=ATG
contig_1	Prodigal_v2.6.3	CDS	401	600	45.1	-	0	ID=PRODIGAL_00002;partial=00;start_type=GTG
"#;
        std::fs::write(&gff_path, content).unwrap();

        let coords = parse_prodigal_gff(&gff_path).unwrap();
        assert_eq!(coords.len(), 2);

        let gene1 = coords.get("PRODIGAL_00001").unwrap();
        assert_eq!(gene1.contig, "contig_1");
        assert_eq!(gene1.start, 101);
        assert_eq!(gene1.end, 300);
        assert_eq!(gene1.strand, Strand::Forward);

        let gene2 = coords.get("PRODIGAL_00002").unwrap();
        assert_eq!(gene2.contig, "contig_1");
        assert_eq!(gene2.start, 401);
        assert_eq!(gene2.end, 600);
        assert_eq!(gene2.strand, Strand::Reverse);
    }

    #[cfg(not(feature = "prodigal"))]
    #[test]
    fn test_orphos_not_available_without_feature() {
        assert!(!OrphosRunner::is_available());
        assert!(!OrphosRunner::is_installed());
    }
}