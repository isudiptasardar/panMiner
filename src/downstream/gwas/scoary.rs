//! Scoary2 gene-trait association testing runner.
//!
//! Scoary2 is a gene-trait association tool that performs Fisher's exact test
//! for each gene cluster against user-supplied phenotypes.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};
use crate::error::{Error, Result};

/// Scoary2Runner performs gene-trait association testing via Scoary2.
pub struct Scoary2Runner {
    scoary2_path: PathBuf,
    phenotypes_file: Option<PathBuf>,
    output_dir: Option<PathBuf>,
    threads: usize,
}

impl Scoary2Runner {
    pub fn new() -> Self {
        Self {
            scoary2_path: PathBuf::from("scoary2"),
            phenotypes_file: None,
            output_dir: None,
            threads: 1,
        }
    }

    pub fn with_phenotypes(mut self, path: PathBuf) -> Self {
        self.phenotypes_file = Some(path);
        self
    }

    pub fn with_output_dir(mut self, path: PathBuf) -> Self {
        self.output_dir = Some(path);
        self
    }

    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    pub fn detect() -> Option<Self> {
        which::which("scoary2").ok().map(|path| Self {
            scoary2_path: path,
            phenotypes_file: None,
            output_dir: None,
            threads: 1,
        })
    }

    fn run_analysis(&self, output_dir: &Path) -> Result<Scoary2Result> {
        let phenotypes_path = self.phenotypes_file.as_ref().ok_or_else(|| {
            Error::Config("Scoary2 requires a phenotypes file. Use `with_phenotypes()` to set it.".to_string())
        })?;

        if !phenotypes_path.exists() {
            return Err(Error::Config(format!(
                "Phenotypes file not found: {}",
                phenotypes_path.display()
            )));
        }

        let pa_csv_path = output_dir.join("gene_presence_absence.csv");
        if !pa_csv_path.exists() {
            return Err(Error::Config(format!(
                "Gene presence/absence CSV not found at {}. Is this a valid PanMiner output directory?",
                pa_csv_path.display()
            )));
        }

        let (gene_ids, genome_ids, presence_matrix) = self.parse_pa_csv(&pa_csv_path)?;

        if gene_ids.is_empty() {
            return Err(Error::Config("No genes found in presence/absence CSV.".to_string()));
        }
        if genome_ids.is_empty() {
            return Err(Error::Config("No genomes found in presence/absence CSV.".to_string()));
        }

        let temp_dir = std::env::temp_dir().join("panminer_scoary2");
        std::fs::create_dir_all(&temp_dir)?;

        let scoary_pa_path = temp_dir.join("presence_absence.csv");
        self.write_scoary_pa_csv(&scoary_pa_path, &gene_ids, &genome_ids, &presence_matrix)?;

        let phenotypes_content = std::fs::read_to_string(phenotypes_path)?;
        self.validate_phenotypes(&phenotypes_content, &genome_ids)?;

        let downstream_dir = output_dir.join("downstream");
        std::fs::create_dir_all(&downstream_dir)?;
        let results_dir = temp_dir.join("results");

        let mut cmd = Command::new(&self.scoary2_path);
        cmd.arg("-t").arg(phenotypes_path)
           .arg("-p").arg(&scoary_pa_path)
           .arg("-o").arg(&results_dir)
           .arg("-e");

        if self.threads > 1 {
            cmd.arg("-T").arg(self.threads.to_string());
        }

        let output = cmd.output().map_err(|e| {
            Error::Config(format!("Failed to execute scoary2: {}. Is scoary2 installed? (pip install scoary2)", e))
        })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::Config(format!(
                "Scoary2 failed with exit code {:?}: {}",
                output.status.code(),
                stderr
            )));
        }

        let results_csv = self.find_scoary_results(&results_dir)?;
        let associations = self.parse_scoary_output(&results_csv)?;
        let output_path = downstream_dir.join("scoary_results.csv");
        self.write_results_csv(&output_path, &associations)?;
        let trait_name = self.extract_trait_name(phenotypes_path)?;

        Ok(Scoary2Result { associations, trait_name })
    }

    fn parse_pa_csv(&self, path: &Path) -> Result<(Vec<String>, Vec<String>, Vec<Vec<u8>>)> {
        let mut rdr = csv::Reader::from_path(path)?;
        let headers = rdr.headers()?.clone();

        if headers.len() < 4 {
            return Err(Error::Config(format!(
                "Invalid P/A CSV format: expected at least 4 columns, got {}",
                headers.len()
            )));
        }

        let genome_ids: Vec<String> = headers.iter().skip(3).map(|s| s.to_string()).collect();
        let mut gene_ids = Vec::new();
        let mut presence_matrix = Vec::new();

        for result in rdr.records() {
            let record = result?;
            if record.len() < 4 {
                continue;
            }
            let gene_id = record.get(0).unwrap_or("").to_string();
            gene_ids.push(gene_id);

            let presence: Vec<u8> = record.iter().skip(3).map(|cell| {
                let val = cell.trim();
                if val.is_empty() || val == "0" { 0 } else { 1 }
            }).collect();
            presence_matrix.push(presence);
        }

        Ok((gene_ids, genome_ids, presence_matrix))
    }

    fn write_scoary_pa_csv(&self, path: &Path, gene_ids: &[String], genome_ids: &[String], presence_matrix: &[Vec<u8>]) -> Result<()> {
        let mut wtr = csv::Writer::from_path(path)?;
        let mut header = vec!["gene".to_string()];
        header.extend(genome_ids.iter().cloned());
        wtr.write_record(&header)?;

        for (gene_id, presence) in gene_ids.iter().zip(presence_matrix.iter()) {
            let mut row = vec![gene_id.clone()];
            for &p in presence {
                row.push(if p == 0 { "0".to_string() } else { "1".to_string() });
            }
            wtr.write_record(&row)?;
        }
        wtr.flush()?;
        Ok(())
    }

    fn validate_phenotypes(&self, content: &str, _expected_genomes: &[String]) -> Result<()> {
        let first_line = content.lines().next().unwrap_or_default();
        let fields: Vec<&str> = if first_line.contains('\t') {
            first_line.split('\t').collect()
        } else if first_line.contains(',') {
            first_line.split(',').collect()
        } else {
            return Err(Error::Config(format!("Phenotypes file must be TSV or CSV format, got: {}", first_line)));
        };

        if fields.len() < 2 {
            return Err(Error::Config("Phenotypes file must have at least 2 columns (genome, phenotype)".to_string()));
        }
        Ok(())
    }

    fn find_scoary_results(&self, results_dir: &Path) -> Result<PathBuf> {
        let entries = std::fs::read_dir(results_dir)
            .map_err(|e| Error::Config(format!("Failed to read scoary2 results directory: {}", e)))?;

        for entry in entries.flatten() {
            let path = entry.path();
            if let Some(name) = path.file_name() {
                let name_str = name.to_string_lossy();
                if name_str.ends_with("_results.csv") || name_str == "results.csv" {
                    return Ok(path);
                }
            }
        }

        let entries = std::fs::read_dir(results_dir)?;
        for entry in entries.flatten() {
            if entry.path().extension().map(|e| e == "csv").unwrap_or(false) {
                return Ok(entry.path());
            }
        }

        Err(Error::Config(format!("No results CSV found in scoary2 output directory: {}", results_dir.display())))
    }

    fn parse_scoary_output(&self, path: &Path) -> Result<Vec<ScoaryAssociation>> {
        let mut rdr = csv::Reader::from_path(path)?;
        let headers = rdr.headers()?.clone();

        let gene_idx = self.find_column(&headers, &["gene", "Gene", "gene_id", "cluster_id"])?;
        let trait_idx = self.find_column(&headers, &["trait", "Trait", "phenotype"])?;
        let pval_idx = self.find_column(&headers, &["p_value", "pvalue", "p_value_corrected"])?;
        let fdr_idx = self.find_column(&headers, &["FDR", "fdr", "padj", "q_value"])?;
        let effect_idx = self.find_column(&headers, &["effect_size", "effect", "odds_ratio"])?;
        let n_present_idx = self.find_column(&headers, &["n_present", "present", "cases"])?;
        let n_absent_idx = self.find_column(&headers, &["n_absent", "absent", "controls"])?;

        let mut associations = Vec::new();
        for result in rdr.records() {
            let record = result?;
            associations.push(ScoaryAssociation {
                gene: record.get(gene_idx).unwrap_or("").to_string(),
                trait_name: record.get(trait_idx).unwrap_or("").to_string(),
                p_value: record.get(pval_idx).and_then(|s| s.parse::<f64>().ok()).unwrap_or(1.0),
                fdr: record.get(fdr_idx).and_then(|s| s.parse::<f64>().ok()).unwrap_or(1.0),
                effect_size: record.get(effect_idx).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0),
                n_present: record.get(n_present_idx).and_then(|s| s.parse::<usize>().ok()).unwrap_or(0),
                n_absent: record.get(n_absent_idx).and_then(|s| s.parse::<usize>().ok()).unwrap_or(0),
            });
        }

        associations.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap_or(std::cmp::Ordering::Equal));
        Ok(associations)
    }

    fn find_column(&self, headers: &csv::StringRecord, possible_names: &[&str]) -> Result<usize> {
        for name in possible_names {
            for (i, header) in headers.iter().enumerate() {
                if header == *name {
                    return Ok(i);
                }
            }
        }
        Err(Error::Config(format!(
            "Could not find any of columns {:?} in CSV header",
            possible_names
        )))
    }

    fn extract_trait_name(&self, phenotypes_path: &Path) -> Result<String> {
        let content = std::fs::read_to_string(phenotypes_path)?;
        let first_line = content.lines().next().unwrap_or_default();
        let fields: Vec<&str> = if first_line.contains('\t') {
            first_line.split('\t').collect()
        } else {
            first_line.split(',').collect()
        };

        if fields.len() >= 2 {
            let first_col = fields[0].trim().to_lowercase();
            if first_col == "genome" || first_col == "sample" || first_col == "id" {
                return Ok(fields[1].trim().to_string());
            }
        }

        Ok(phenotypes_path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_else(|| "unknown_trait".to_string()))
    }

    fn write_results_csv(&self, path: &Path, associations: &[ScoaryAssociation]) -> Result<()> {
        let mut wtr = csv::Writer::from_path(path)?;
        wtr.write_record(&["gene", "trait", "p_value", "FDR", "effect_size", "n_present", "n_absent"])?;
        for assoc in associations {
            wtr.write_record(&[
                &assoc.gene, &assoc.trait_name, &assoc.p_value.to_string(),
                &assoc.fdr.to_string(), &assoc.effect_size.to_string(),
                &assoc.n_present.to_string(), &assoc.n_absent.to_string(),
            ])?;
        }
        wtr.flush()?;
        Ok(())
    }
}

impl Default for Scoary2Runner {
    fn default() -> Self { Self::new() }
}

impl DownstreamRunner for Scoary2Runner {
    type Output = Scoary2Result;

    fn run(&self, output_dir: &Path) -> Result<Self::Output> {
        self.run_analysis(output_dir)
    }

    fn name(&self) -> &str { "Scoary2" }

    fn is_available(&self) -> bool { Self::detect().is_some() }

    fn required_inputs(&self) -> Vec<DownstreamInput> {
        vec![DownstreamInput::PresenceAbsenceCsv, DownstreamInput::PhenotypesFile]
    }
}

pub struct Scoary2Result {
    pub associations: Vec<ScoaryAssociation>,
    pub trait_name: String,
}

impl DownstreamResult for Scoary2Result {
    fn write_to(&self, dir: &Path) -> Result<()> {
        std::fs::create_dir_all(dir)?;
        let csv_path = dir.join("scoary_results.csv");
        let mut wtr = csv::Writer::from_path(&csv_path)?;
        wtr.write_record(&["gene", "trait", "p_value", "FDR", "effect_size", "n_present", "n_absent"])?;
        for assoc in &self.associations {
            wtr.write_record(&[
                &assoc.gene, &assoc.trait_name, &assoc.p_value.to_string(),
                &assoc.fdr.to_string(), &assoc.effect_size.to_string(),
                &assoc.n_present.to_string(), &assoc.n_absent.to_string(),
            ])?;
        }
        wtr.flush()?;
        Ok(())
    }

    fn summary(&self) -> String {
        let n_associations = self.associations.len();
        let significant = self.associations.iter().filter(|a| a.fdr < 0.05).count();
        format!("Scoary2: trait={}, total_associations={}, significant_fdr_0.05={}",
            self.trait_name, n_associations, significant)
    }
}

#[derive(Debug, Clone)]
pub struct ScoaryAssociation {
    pub gene: String,
    pub trait_name: String,
    pub p_value: f64,
    pub fdr: f64,
    pub effect_size: f64,
    pub n_present: usize,
    pub n_absent: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    fn make_test_pa_csv(path: &Path) -> Result<()> {
        let mut wtr = csv::Writer::from_path(path)?;
        wtr.write_record(&["gene", "annotation", "cluster_id", "genome1", "genome2", "genome3", "genome4"])?;
        wtr.write_record(&["gene_001", "hypothetical protein", "cluster_001", "1", "1", "0", "0"])?;
        wtr.write_record(&["gene_002", "transposase", "cluster_002", "1", "0", "1", "0"])?;
        wtr.write_record(&["gene_003", "ABC transporter", "cluster_003", "0", "1", "1", "1"])?;
        wtr.flush()?;
        Ok(())
    }

    fn make_test_phenotypes(path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(path)?;
        writeln!(file, "genome\tphenotype")?;
        writeln!(file, "genome1\t1")?;
        writeln!(file, "genome2\t0")?;
        writeln!(file, "genome3\t1")?;
        writeln!(file, "genome4\t0")?;
        Ok(())
    }

    #[test]
    fn test_parse_pa_csv() {
        let temp_dir = TempDir::new().unwrap();
        let pa_path = temp_dir.path().join("pa.csv");
        make_test_pa_csv(&pa_path).unwrap();
        let runner = Scoary2Runner::new();
        let (gene_ids, genome_ids, presence) = runner.parse_pa_csv(&pa_path).unwrap();
        assert_eq!(gene_ids.len(), 3);
        assert_eq!(genome_ids.len(), 4);
        assert_eq!(presence.len(), 3);
        assert_eq!(presence[0], vec![1, 1, 0, 0]);
        assert_eq!(presence[1], vec![1, 0, 1, 0]);
        assert_eq!(presence[2], vec![0, 1, 1, 1]);
    }

    #[test]
    fn test_write_scoary_pa_csv() {
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().join("scoary_pa.csv");
        let runner = Scoary2Runner::new();
        let gene_ids = vec!["gene_001".to_string(), "gene_002".to_string()];
        let genome_ids = vec!["genome1".to_string(), "genome2".to_string(), "genome3".to_string()];
        let presence = vec![vec![1, 1, 0], vec![1, 0, 1]];
        runner.write_scoary_pa_csv(&output_path, &gene_ids, &genome_ids, &presence).unwrap();
        let content = std::fs::read_to_string(&output_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines[0], "gene,genome1,genome2,genome3");
        assert_eq!(lines[1], "gene_001,1,1,0");
        assert_eq!(lines[2], "gene_002,1,0,1");
    }

    #[test]
    fn test_extract_trait_name_from_header() {
        let temp_dir = TempDir::new().unwrap();
        let pheno_path = temp_dir.path().join("phenotypes.txt");
        make_test_phenotypes(&pheno_path).unwrap();
        let runner = Scoary2Runner::new();
        let trait_name = runner.extract_trait_name(&pheno_path).unwrap();
        assert_eq!(trait_name, "phenotype");
    }

    #[test]
    fn test_scoary2_result_summary() {
        let result = Scoary2Result {
            associations: vec![
                ScoaryAssociation { gene: "gene_001".to_string(), trait_name: "resistance".to_string(), p_value: 0.001, fdr: 0.01, effect_size: 2.5, n_present: 10, n_absent: 5 },
                ScoaryAssociation { gene: "gene_002".to_string(), trait_name: "resistance".to_string(), p_value: 0.05, fdr: 0.08, effect_size: 1.2, n_present: 8, n_absent: 7 },
            ],
            trait_name: "resistance".to_string(),
        };
        let summary = result.summary();
        assert!(summary.contains("Scoary2"));
        assert!(summary.contains("resistance"));
        assert!(summary.contains("total_associations=2"));
        assert!(summary.contains("significant_fdr_0.05=1"));
    }

    #[test]
    fn test_scoary2_result_write_to() {
        let temp_dir = TempDir::new().unwrap();
        let result_dir = temp_dir.path().join("results");
        let result = Scoary2Result {
            associations: vec![ScoaryAssociation { gene: "gene_001".to_string(), trait_name: "resistance".to_string(), p_value: 0.001, fdr: 0.01, effect_size: 2.5, n_present: 10, n_absent: 5 }],
            trait_name: "resistance".to_string(),
        };
        result.write_to(&result_dir).unwrap();
        let csv_path = result_dir.join("scoary_results.csv");
        assert!(csv_path.exists());
        let content = std::fs::read_to_string(&csv_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines[0], "gene,trait,p_value,FDR,effect_size,n_present,n_absent");
        assert!(lines[1].contains("gene_001"));
        assert!(lines[1].contains("0.001"));
    }
}
