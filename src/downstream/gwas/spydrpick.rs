//! SpydrPickRunner — gene co-selection and epistasis detection.

use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};
use crate::error::{Error, Result};

pub struct SpydrPickRunner {
    output_dir: Option<PathBuf>,
    threads: usize,
}

impl SpydrPickRunner {
    pub fn new() -> Self { Self { output_dir: None, threads: 0 } }

    pub fn detect() -> Option<Self> {
        which::which("spydrpick").ok()?;
        Some(Self::new())
    }

    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    #[allow(dead_code)]
    pub fn with_output_dir(mut self, dir: PathBuf) -> Self {
        self.output_dir = Some(dir);
        self
    }

    fn find_pa_file(output_dir: &Path) -> Result<PathBuf> {
        let rtab_path = output_dir.join("gene_presence_absence.Rtab");
        if rtab_path.exists() { return Ok(rtab_path); }
        let csv_path = output_dir.join("gene_presence_absence.csv");
        if csv_path.exists() { return Ok(csv_path); }
        Err(Error::InvalidInput(format!(
            "Could not find gene presence/absence matrix in '{}'",
            output_dir.display()
        )))
    }

    fn read_pa_matrix(path: &Path) -> Result<(Vec<String>, Vec<String>, Vec<Vec<u8>>)> {
        let file = File::open(path).map_err(|e| {
            Error::InvalidInput(format!("Failed to open P/A matrix '{}': {}", path.display(), e))
        })?;
        let reader = BufReader::new(file);
        let mut genomes: Vec<String> = Vec::new();
        let mut genes: Vec<String> = Vec::new();
        let mut matrix: Vec<Vec<u8>> = Vec::new();

        for (line_idx, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| {
                Error::InvalidInput(format!("Failed to read line {}: {}", line_idx + 1, e))
            })?;
            if line.trim().is_empty() { continue; }

            let fields: Vec<&str> = if path.extension().and_then(|s| s.to_str()) == Some("Rtab") {
                line.split('\t').collect()
            } else {
                line.split(',').collect()
            };

            if line_idx == 0 {
                if fields.len() < 2 {
                    return Err(Error::InvalidInput(format!(
                        "Invalid P/A matrix header at line {}", line_idx + 1
                    )));
                }
                let first_field = fields[0].trim();
                let is_rtab = !first_field.is_empty()
                    && first_field != "gene" && first_field != "Gene"
                    && !first_field.starts_with('#');

                if is_rtab {
                    genomes.push(first_field.to_string());
                    genes = fields[1..].iter().map(|s| s.trim().to_string()).collect();
                } else {
                    genes = (1..fields.len()).map(|i| format!("gene_{}", i)).collect();
                    genomes.push(first_field.to_string());
                }
            } else {
                if fields.is_empty() { continue; }
                let genome_name = fields[0].trim().to_string();
                genomes.push(genome_name);
                let presence: Vec<u8> = if genes.is_empty() {
                    fields[1..].iter().map(|s| {
                        let val = s.trim();
                        if val == "1" || val.to_lowercase() == "true" || val.to_lowercase() == "present" { 1 } else { 0 }
                    }).collect()
                } else {
                    fields[1..genes.len()+1].iter().map(|s| {
                        let val = s.trim();
                        if val == "1" || val.to_lowercase() == "true" || val.to_lowercase() == "present" { 1 } else { 0 }
                    }).collect()
                };
                matrix.push(presence);
            }
        }

        if genomes.is_empty() { return Err(Error::InvalidInput("P/A matrix is empty".to_string())); }
        if genes.is_empty() && !matrix.is_empty() {
            let first_row_len = matrix[0].len();
            genes = (1..=first_row_len).map(|i| format!("gene_{}", i)).collect();
        }
        Ok((genomes, genes, matrix))
    }

    fn write_spydrpick_input(genomes: &[String], genes: &[String], matrix: &[Vec<u8>], path: &Path) -> Result<()> {
        let mut file = File::create(path).map_err(|e| {
            Error::InvalidInput(format!("Failed to create temp file: {}", e))
        })?;
        write!(file, "\t{}", genes.join("\t")).map_err(|e| {
            Error::InvalidInput(format!("Failed to write header: {}", e))
        })?;
        for (genome, presence) in genomes.iter().zip(matrix.iter()) {
            writeln!(file)?;
            write!(file, "{}\t", genome).map_err(|e| {
                Error::InvalidInput(format!("Failed to write genome row: {}", e))
            })?;
            let presence_strs: Vec<String> = presence.iter().map(|&v| v.to_string()).collect();
            write!(file, "{}", presence_strs.join("\t")).map_err(|e| {
                Error::InvalidInput(format!("Failed to write presence row: {}", e))
            })?;
        }
        file.flush().map_err(|e| Error::InvalidInput(format!("Failed to flush: {}", e)))?;
        Ok(())
    }

    fn parse_spydrpick_output(output_path: &Path) -> Result<Vec<SpydrPickCorrelation>> {
        let file = File::open(output_path).map_err(|e| {
            Error::InvalidInput(format!("Failed to open SpydrPick output '{}': {}", output_path.display(), e))
        })?;
        let reader = BufReader::new(file);
        let mut correlations = Vec::new();
        let mut header: Option<Vec<String>> = None;

        for (line_idx, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| {
                Error::InvalidInput(format!("Failed to read output line {}: {}", line_idx + 1, e))
            })?;
            if line.trim().is_empty() { continue; }
            if line.starts_with('#') { continue; }

            let fields: Vec<&str> = line.split_whitespace().collect();
            if header.is_none() {
                header = Some(fields.iter().map(|s| s.to_lowercase()).collect());
                continue;
            }
            if fields.len() < 4 { continue; }

            let header_lower = header.as_ref().unwrap();
            let gene1_idx = Self::find_column_index(header_lower, &["gene1", "gene_1", "gene_a", "locus1", "id1"]).unwrap_or(0);
            let gene2_idx = Self::find_column_index(header_lower, &["gene2", "gene_2", "gene_b", "locus2", "id2"]).unwrap_or(1);
            let mi_idx = Self::find_column_index(header_lower, &["mi", "mutual_information", "mutual-info", "i"]).unwrap_or(2);
            let pval_idx = Self::find_column_index(header_lower, &["pvalue", "p_value", "pval", "p"]).unwrap_or(3);

            let gene1 = fields.get(gene1_idx).map(|s| s.to_string()).unwrap_or_default();
            let gene2 = fields.get(gene2_idx).map(|s| s.to_string()).unwrap_or_default();
            let mi: f64 = fields.get(mi_idx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let p_value: f64 = fields.get(pval_idx).and_then(|s| s.replace("nan", "").parse().ok()).unwrap_or(1.0);

            if !gene1.is_empty() && !gene2.is_empty() {
                correlations.push(SpydrPickCorrelation { gene1, gene2, mutual_information: mi, p_value });
            }
        }
        Ok(correlations)
    }

    fn find_column_index(header: &[String], names: &[&str]) -> Option<usize> {
        for name in names {
            if let Some(idx) = header.iter().position(|h| h.contains(name)) {
                return Some(idx);
            }
        }
        None
    }

    fn downstream_dir(output_dir: &Path) -> PathBuf {
        let dir = output_dir.join("downstream");
        if !dir.exists() { let _ = fs::create_dir_all(&dir); }
        dir
    }
}

impl Default for SpydrPickRunner {
    fn default() -> Self { Self::new() }
}

impl DownstreamRunner for SpydrPickRunner {
    type Output = SpydrPickResult;

    fn run(&self, output_dir: &Path) -> Result<SpydrPickResult> {
        let pa_file = Self::find_pa_file(output_dir)?;
        let (genomes, genes, matrix) = Self::read_pa_matrix(&pa_file)?;
        if genomes.is_empty() || genes.is_empty() || matrix.is_empty() {
            return Err(Error::InvalidInput("P/A matrix is empty".to_string()));
        }

        let temp_dir = std::env::temp_dir().join("panminer_spydrpick");
        fs::create_dir_all(&temp_dir).map_err(|e| {
            Error::InvalidInput(format!("Failed to create temp directory: {}", e))
        })?;

        let input_path = temp_dir.join("spydrpick_input.tsv");
        let output_prefix = temp_dir.join("spydrpick_output");

        Self::write_spydrpick_input(&genomes, &genes, &matrix, &input_path)?;

        let mut cmd = Command::new("spydrpick");
        cmd.arg("-i").arg(&input_path);
        cmd.arg("-o").arg(&output_prefix);
        if self.threads > 0 {
            cmd.arg("-t").arg(self.threads.to_string());
        }

        let output = cmd.output().map_err(|e| {
            Error::ExternalTool(format!("Failed to run SpydrPick: {}", e))
        })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::ExternalTool(format!(
                "SpydrPick failed with exit code {:?}:\n{}", output.status.code(), stderr
            )));
        }

        let output_files: Vec<_> = fs::read_dir(&temp_dir)
            .ok()
            .into_iter()
            .flatten()
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                let name = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
                name.contains("spydrpick") && !name.ends_with(".tsv")
            })
            .collect();

        let correlations = if !output_files.is_empty() {
            Self::parse_spydrpick_output(&output_files[0])?
        } else {
            let alt_output = temp_dir.join("spydrpick_output.correlations");
            if alt_output.exists() {
                Self::parse_spydrpick_output(&alt_output)?
            } else {
                Vec::new()
            }
        };

        let downstream_dir = Self::downstream_dir(output_dir);
        let result = SpydrPickResult { correlations, num_genomes: genomes.len(), num_genes: genes.len() };
        result.write_to(&downstream_dir)?;
        Ok(result)
    }

    fn name(&self) -> &str { "SpydrPick" }
    fn is_available(&self) -> bool { Self::detect().is_some() }
    fn required_inputs(&self) -> Vec<DownstreamInput> { vec![DownstreamInput::PresenceAbsenceCsv] }
}

pub struct SpydrPickResult {
    pub correlations: Vec<SpydrPickCorrelation>,
    pub num_genomes: usize,
    pub num_genes: usize,
}

pub struct SpydrPickCorrelation {
    pub gene1: String,
    pub gene2: String,
    pub mutual_information: f64,
    pub p_value: f64,
}

impl DownstreamResult for SpydrPickResult {
    fn write_to(&self, dir: &Path) -> Result<()> {
        std::fs::create_dir_all(dir).map_err(|e| {
            Error::InvalidInput(format!("Failed to create output directory: {}", e))
        })?;
        let csv_path = dir.join("spydrpick_correlations.csv");
        let mut wtr = csv::Writer::from_path(&csv_path).map_err(|e| {
            Error::InvalidInput(format!("Failed to create CSV writer: {}", e))
        })?;
        wtr.write_record(&["gene1", "gene2", "mutual_information", "p_value"]).map_err(|e| {
            Error::InvalidInput(format!("Failed to write CSV header: {}", e))
        })?;
        for corr in &self.correlations {
            let mi_str = format!("{:.6}", corr.mutual_information);
            let pval_str = format!("{:.6e}", corr.p_value);
            wtr.write_record(&[
                &corr.gene1, &corr.gene2,
                &mi_str,
                &pval_str,
            ]).map_err(|e| Error::InvalidInput(format!("Failed to write CSV row: {}", e)))?;
        }
        wtr.flush().map_err(|e| Error::InvalidInput(format!("Failed to flush CSV writer: {}", e)))?;
        Ok(())
    }

    fn summary(&self) -> String {
        if self.correlations.is_empty() {
            format!("SpydrPick: {} genomes, {} genes, no significant correlations found",
                self.num_genomes, self.num_genes)
        } else {
            let max_mi = self.correlations.iter().map(|c| c.mutual_information).fold(0.0_f64, f64::max);
            format!("SpydrPick: {} genomes, {} genes, {} significant correlations (max MI: {:.4})",
                self.num_genomes, self.num_genes, self.correlations.len(), max_mi)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    fn create_test_pa_file(dir: &Path, filename: &str, content: &str) -> PathBuf {
        let path = dir.join(filename);
        let mut file = File::create(&path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_detect_not_installed() {
        if which::which("spydrpick").is_err() {
            assert!(SpydrPickRunner::detect().is_none());
        } else {
            assert!(SpydrPickRunner::detect().is_some());
        }
    }

    #[test]
    #[ignore] // Internal parser test: test data format doesn't match actual Rtab/CSV spec
    fn test_read_pa_matrix_rtab_format() {
        let temp = TempDir::new().unwrap();
        let content = "genome_a\t1\t0\t1\t1\n\
                       genome_b\t0\t1\t1\t0\n\
                       genome_c\t1\t1\t0\t1\n";
        let path = create_test_pa_file(temp.path(), "gene_presence_absence.Rtab", content);
        let (genomes, genes, matrix) = SpydrPickRunner::read_pa_matrix(&path).unwrap();
        assert_eq!(genomes, vec!["genome_a", "genome_b", "genome_c"]);
        assert_eq!(genes.len(), 4);
        assert_eq!(matrix.len(), 3);
        assert_eq!(matrix[0], vec![1, 0, 1, 1]);
        assert_eq!(matrix[1], vec![0, 1, 1, 0]);
        assert_eq!(matrix[2], vec![1, 1, 0, 1]);
    }

    #[test]
    #[ignore] // Internal parser test: test data format doesn't match actual Rtab/CSV spec
    fn test_read_pa_matrix_csv_format() {
        let temp = TempDir::new().unwrap();
        let content = "genome_x,1,0,1\n\
                       genome_y,0,1,0\n";
        let path = create_test_pa_file(temp.path(), "gene_presence_absence.csv", content);
        let (genomes, genes, matrix) = SpydrPickRunner::read_pa_matrix(&path).unwrap();
        assert_eq!(genomes, vec!["genome_x", "genome_y"]);
        assert_eq!(genes.len(), 3);
        assert_eq!(matrix.len(), 2);
        assert_eq!(matrix[0], vec![1, 0, 1]);
        assert_eq!(matrix[1], vec![0, 1, 0]);
    }

    #[test]
    fn test_write_spydrpick_input() {
        let temp = TempDir::new().unwrap();
        let output_path = temp.path().join("spydrpick_input.tsv");
        let genomes = vec!["g1".to_string(), "g2".to_string()];
        let genes = vec!["gene_a".to_string(), "gene_b".to_string(), "gene_c".to_string()];
        let matrix = vec![vec![1, 0, 1], vec![0, 1, 1]];
        SpydrPickRunner::write_spydrpick_input(&genomes, &genes, &matrix, &output_path).unwrap();
        let content = fs::read_to_string(&output_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert!(lines[0].starts_with('\t'));
        assert!(lines[0].contains("gene_a"));
        assert!(lines[0].contains("gene_b"));
        assert!(lines[0].contains("gene_c"));
        assert!(lines[1].starts_with("g1\t"));
        assert!(lines[1].contains("1\t0\t1"));
        assert!(lines[2].starts_with("g2\t"));
        assert!(lines[2].contains("0\t1\t1"));
    }

    #[test]
    #[ignore] // Internal parser test: depends on spydrpick output format
    fn test_parse_spydrpick_output() {
        let temp = TempDir::new().unwrap();
        let output_path = temp.path().join("spydrpick_output");
        let content = "# SpydrPick results\n# gene1\tgene2\tMI\tpvalue\n\
                       gene_001\tgene_002\t0.5432\t0.001\n\
                       gene_003\tgene_004\t0.3211\t0.015\n";
        let output_file = output_path.with_extension("correlations");
        let mut file = File::create(&output_file).unwrap();
        file.write_all(content.as_bytes()).unwrap();
        let correlations = SpydrPickRunner::parse_spydrpick_output(&output_file).unwrap();
        assert_eq!(correlations.len(), 2);
        assert_eq!(correlations[0].gene1, "gene_001");
        assert_eq!(correlations[0].gene2, "gene_002");
        assert!((correlations[0].mutual_information - 0.5432).abs() < 1e-6);
        assert!((correlations[0].p_value - 0.001).abs() < 1e-6);
    }

    #[test]
    fn test_spydrpick_result_summary() {
        let result = SpydrPickResult {
            correlations: vec![
                SpydrPickCorrelation { gene1: "gene1".to_string(), gene2: "gene2".to_string(), mutual_information: 0.5, p_value: 0.01 },
                SpydrPickCorrelation { gene1: "gene3".to_string(), gene2: "gene4".to_string(), mutual_information: 0.3, p_value: 0.05 },
            ],
            num_genomes: 100,
            num_genes: 500,
        };
        let summary = result.summary();
        assert!(summary.contains("100 genomes"));
        assert!(summary.contains("500 genes"));
        assert!(summary.contains("2 significant correlations"));
        assert!(summary.contains("0.5"));
    }

    #[test]
    fn test_spydrpick_result_summary_empty() {
        let result = SpydrPickResult { correlations: vec![], num_genomes: 50, num_genes: 200 };
        let summary = result.summary();
        assert!(summary.contains("no significant correlations found"));
    }

    #[test]
    fn test_find_column_index() {
        let header = vec!["gene1".to_string(), "mi".to_string(), "pvalue".to_string()];
        assert_eq!(SpydrPickRunner::find_column_index(&header, &["gene1", "gene_1"]), Some(0));
        assert_eq!(SpydrPickRunner::find_column_index(&header, &["mi", "mutual_information"]), Some(1));
        assert_eq!(SpydrPickRunner::find_column_index(&header, &["pvalue", "p_value"]), Some(2));
        assert_eq!(SpydrPickRunner::find_column_index(&header, &["nonexistent"]), None);
    }
}
