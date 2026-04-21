//! AmrFinderRunner — AMR gene detection via AMRFinderPlus.
//!
//! AMRFinderPlus provides NCBI-curated AMR detection with hierarchical
//! evidence levels: EXACT > ALLELE > BLAST > HMM > PARTIAL > POINT > INTERNAL_STOP.

use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};
use crate::error::{Error, Result};

/// AMRFinderPlus detection runner.
pub struct AmrFinderRunner {
    database_path: Option<PathBuf>,
    organism: Option<String>,
    threads: usize,
}

impl AmrFinderRunner {
    pub fn new() -> Self {
        Self {
            database_path: None,
            organism: None,
            threads: 1,
        }
    }

    pub fn with_database(mut self, path: PathBuf) -> Self {
        self.database_path = Some(path);
        self
    }

    pub fn with_organism(mut self, org: String) -> Self {
        self.organism = Some(org);
        self
    }

    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    pub fn detect() -> Option<Self> {
        which::which("amrfinder")
            .ok()
            .or_else(|| which::which("amrfinder_plus").ok())
            .map(|_| Self::new())
    }

    fn amrfinder_path() -> Option<PathBuf> {
        which::which("amrfinder").ok().or_else(|| which::which("amrfinder_plus").ok())
    }

    fn run_internal(&self, output_dir: &Path) -> Result<AmrFinderResult> {
        let amrfinder_path = Self::amrfinder_path()
            .ok_or_else(|| Error::ExternalTool("AMRFinderPlus not found".to_string()))?;

        let temp_dir = tempfile::TempDir::new().map_err(|e| Error::Io(e))?;
        let temp_input_dir = temp_dir.path();
        let temp_output_dir = temp_dir.path();

        let protein_fasta_path = output_dir.join("combined_protein_CDS.fasta");
        if !protein_fasta_path.exists() {
            return Err(Error::InvalidInput(format!(
                "Protein FASTA not found at {:?}",
                protein_fasta_path
            )));
        }

        let input_fasta = temp_input_dir.join("input_proteins.fasta");
        fs::copy(&protein_fasta_path, &input_fasta).map_err(|e| Error::Io(e))?;

        let gene_data_path = output_dir.join("gene_data.csv");
        let gene_data_gff = if gene_data_path.exists() {
            let gff_path = temp_input_dir.join("gene_data.gff");
            convert_csv_to_gff(&gene_data_path, &gff_path)?;
            Some(gff_path)
        } else {
            None
        };

        let output_prefix = temp_output_dir.join("amr_output");
        let mut cmd = Command::new(&amrfinder_path);
        cmd.arg("-i").arg(&input_fasta);
        cmd.arg("-o").arg(&output_prefix);
        cmd.arg("--plus");

        if let Some(ref db_path) = self.database_path {
            cmd.arg("-d").arg(db_path);
        }

        if let Some(ref organism) = self.organism {
            cmd.arg("--organism").arg(organism);
        }

        if self.threads > 1 {
            cmd.arg("-t").arg(self.threads.to_string());
        }

        if let Some(ref gff_path) = gene_data_gff {
            cmd.arg("--gff").arg(gff_path);
        }

        let output = cmd.output().map_err(|e| Error::ExternalTool(format!(
            "Failed to run AMRFinderPlus: {}", e
        )))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::ExternalTool(format!("AMRFinderPlus failed: {}", stderr)));
        }

        let amr_tsv_path = output_prefix.with_extension("tsv");
        let genes = if amr_tsv_path.exists() {
            parse_amrfinder_tsv(&amr_tsv_path)?
        } else {
            parse_amrfinder_tsv(&output_prefix)?
        };

        let summary_by_class = build_summary(&genes);

        Ok(AmrFinderResult { genes, summary_by_class })
    }
}

impl Default for AmrFinderRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl DownstreamRunner for AmrFinderRunner {
    type Output = AmrFinderResult;

    fn run(&self, output_dir: &Path) -> Result<Self::Output> {
        self.run_internal(output_dir)
    }

    fn name(&self) -> &str {
        "AMRFinderPlus"
    }

    fn is_available(&self) -> bool {
        Self::detect().is_some()
    }

    fn required_inputs(&self) -> Vec<DownstreamInput> {
        vec![DownstreamInput::ProteinFasta, DownstreamInput::GeneDataCsv]
    }
}

pub struct AmrFinderResult {
    pub genes: Vec<AmrGene>,
    pub summary_by_class: HashMap<String, usize>,
}

#[derive(Debug, Clone)]
pub struct AmrGene {
    pub gene_name: String,
    pub scope: String,
    pub target_type: String,
    pub method: String,
    pub identity: f64,
    pub coverage: f64,
    pub contig: String,
    pub start: usize,
    pub end: usize,
    pub strand: String,
    pub annotation: String,
    pub product: String,
    pub resistance: String,
}

impl DownstreamResult for AmrFinderResult {
    fn write_to(&self, dir: &Path) -> Result<()> {
        fs::create_dir_all(dir)?;
        let tsv_path = dir.join("amr_results.tsv");
        write_amr_tsv(&self.genes, &tsv_path)?;
        let summary_path = dir.join("amr_summary.txt");
        write_amr_summary(&self.genes, &self.summary_by_class, &summary_path)?;
        Ok(())
    }

    fn summary(&self) -> String {
        format!(
            "AMRFinderPlus: {} AMR genes detected across {} drug class(es)",
            self.genes.len(),
            self.summary_by_class.len()
        )
    }
}

fn parse_amrfinder_tsv(path: &Path) -> Result<Vec<AmrGene>> {
    let content = fs::read_to_string(path).map_err(|e| Error::Io(e))?;
    let mut genes = Vec::new();
    let mut lines = content.lines();

    let header_line = match lines.next() {
        Some(h) => h,
        None => return Ok(genes),
    };
    let header_fields: Vec<&str> = header_line.split('\t').collect();

    let col = |names: &[&str]| -> usize {
        for name in names {
            if let Some(idx) = header_fields.iter().position(|h| *h == *name) {
                return idx;
            }
        }
        0
    };

    let gene_name_idx = col(&["Gene symbol", "gene_symbol", "Protein ID", "gene_name"]);
    let scope_idx = col(&["Scope", "scope"]);
    let target_type_idx = col(&["Element type", "element_type", "target_type"]);
    let method_idx = col(&["Method", "method", "Evidence"]);
    let identity_idx = col(&["Sequence identity", "seq_identity", "identity", "Identity"]);
    let coverage_idx = col(&["Coverage", "coverage"]);
    let contig_idx = col(&["Contig id", "contig_id", "contig"]);
    let start_idx = col(&["Start", "start", "Gene start", "gene_start"]);
    let end_idx = col(&["Stop", "stop", "end", "Gene stop", "gene_stop"]);
    let strand_idx = col(&["Strand", "strand"]);
    let annotation_idx = col(&["AA Mutation", "aa_mutation", "annotation"]);
    let product_idx = col(&["Product category", "product_category", "product"]);
    let resistance_idx = col(&["Drug class", "drug_class", "resistance", "Class"]);

    for line in lines {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        genes.push(AmrGene {
            gene_name: fields.get(gene_name_idx).unwrap_or(&"").to_string(),
            scope: fields.get(scope_idx).unwrap_or(&"").to_string(),
            target_type: fields.get(target_type_idx).unwrap_or(&"").to_string(),
            method: fields.get(method_idx).unwrap_or(&"").to_string(),
            identity: fields.get(identity_idx).unwrap_or(&"0.0").parse().unwrap_or(0.0),
            coverage: fields.get(coverage_idx).unwrap_or(&"0.0").parse().unwrap_or(0.0),
            contig: fields.get(contig_idx).unwrap_or(&"").to_string(),
            start: fields.get(start_idx).unwrap_or(&"0").parse().unwrap_or(0),
            end: fields.get(end_idx).unwrap_or(&"0").parse().unwrap_or(0),
            strand: fields.get(strand_idx).unwrap_or(&"").to_string(),
            annotation: fields.get(annotation_idx).unwrap_or(&"").to_string(),
            product: fields.get(product_idx).unwrap_or(&"").to_string(),
            resistance: fields.get(resistance_idx).unwrap_or(&"").to_string(),
        });
    }
    Ok(genes)
}

fn write_amr_tsv(genes: &[AmrGene], path: &Path) -> Result<()> {
    let mut wtr = csv::Writer::from_path(path).map_err(|e| Error::Csv(e))?;
    wtr.write_record(&[
        "gene_name", "scope", "target_type", "method", "identity",
        "coverage", "contig", "start", "end", "strand", "annotation",
        "product", "resistance",
    ])?;
    for gene in genes {
        wtr.write_record(&[
            &gene.gene_name, &gene.scope, &gene.target_type, &gene.method,
            &format!("{:.2}", gene.identity), &format!("{:.2}", gene.coverage),
            &gene.contig, &gene.start.to_string(), &gene.end.to_string(),
            &gene.strand, &gene.annotation, &gene.product, &gene.resistance,
        ])?;
    }
    wtr.flush().map_err(Error::Io)?;
    Ok(())
}

fn build_summary(genes: &[AmrGene]) -> HashMap<String, usize> {
    let mut summary: HashMap<String, usize> = HashMap::new();
    for gene in genes {
        *summary.entry(gene.resistance.clone()).or_insert(0) += 1;
    }
    summary
}

fn write_amr_summary(
    genes: &[AmrGene],
    summary_by_class: &HashMap<String, usize>,
    path: &Path,
) -> Result<()> {
    let mut file = fs::File::create(path).map_err(|e| Error::Io(e))?;
    writeln!(&mut file, "AMRFinderPlus AMR Summary")?;
    writeln!(&mut file, "=========================")?;
    writeln!(&mut file)?;
    writeln!(&mut file, "Total AMR genes detected: {}", genes.len())?;
    writeln!(&mut file)?;

    let mut evidence_counts: HashMap<String, usize> = HashMap::new();
    for gene in genes {
        *evidence_counts.entry(gene.method.clone()).or_insert(0) += 1;
    }
    writeln!(&mut file, "Evidence levels:")?;
    for (method, count) in &evidence_counts {
        writeln!(&mut file, "  {}: {}", method, count)?;
    }
    writeln!(&mut file)?;
    writeln!(&mut file, "Genes by drug class:")?;
    let mut sorted_classes: Vec<_> = summary_by_class.iter().collect();
    sorted_classes.sort_by(|a, b| b.1.cmp(a.1));
    for (class, count) in sorted_classes {
        writeln!(&mut file, "  {}: {}", class, count)?;
    }
    writeln!(&mut file)?;
    writeln!(&mut file, "Detailed results:")?;
    for gene in genes {
        writeln!(
            &mut file,
            "  {} [{}] - {} (method: {}, identity: {:.1}%, coverage: {:.1}%)",
            gene.gene_name, gene.resistance, gene.product, gene.method,
            gene.identity, gene.coverage
        )?;
    }
    Ok(())
}

fn convert_csv_to_gff(csv_path: &Path, gff_path: &Path) -> Result<()> {
    let content = fs::read_to_string(csv_path).map_err(|e| Error::Io(e))?;
    let mut file = fs::File::create(gff_path).map_err(|e| Error::Io(e))?;
    let mut lines = content.lines();

    if let Some(_header) = lines.next() {
        writeln!(&mut file, "##gff-version 3")?;
    }

    for line in lines {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 7 {
            continue;
        }
        let cluster_id = fields.get(0).unwrap_or(&"");
        let annotation = fields.get(1).unwrap_or(&"");
        let product = fields.get(2).unwrap_or(&"");
        let contig = fields.get(3).unwrap_or(&"unknown");
        let start = fields.get(4).unwrap_or(&"0");
        let end = fields.get(5).unwrap_or(&"0");
        let strand = fields.get(6).unwrap_or(&"+");

        writeln!(
            &mut file,
            "{}\tpanminer\tgene\t{}\t{}\t.\t{}\t.\tID={};product={};gene={}",
            contig, start, end, strand, cluster_id,
            product.replace(';', "&semicolon;").replace('=', "&equals;"),
            annotation.replace(';', "&semicolon;").replace('=', "&equals;")
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_summary() {
        let genes = vec![
            AmrGene {
                gene_name: "TEM-1".to_string(), scope: "CORE".to_string(),
                target_type: "GENE".to_string(), method: "EXACT".to_string(),
                identity: 100.0, coverage: 100.0, contig: "contig1".to_string(),
                start: 100, end: 200, strand: "+".to_string(),
                annotation: "TEM-1 beta-lactamase".to_string(),
                product: "beta-lactamase".to_string(), resistance: "beta-lactam".to_string(),
            },
            AmrGene {
                gene_name: "aac(3)-IIa".to_string(), scope: "ACCESSORY".to_string(),
                target_type: "GENE".to_string(), method: "ALLELE".to_string(),
                identity: 99.5, coverage: 100.0, contig: "contig2".to_string(),
                start: 300, end: 400, strand: "+".to_string(),
                annotation: "aac(3)-IIa".to_string(),
                product: "aminoglycoside acetyltransferase".to_string(),
                resistance: "aminoglycoside".to_string(),
            },
            AmrGene {
                gene_name: "TEM-2".to_string(), scope: "CORE".to_string(),
                target_type: "GENE".to_string(), method: "BLAST".to_string(),
                identity: 95.0, coverage: 90.0, contig: "contig3".to_string(),
                start: 500, end: 600, strand: "+".to_string(),
                annotation: "TEM-2 beta-lactamase".to_string(),
                product: "beta-lactamase".to_string(), resistance: "beta-lactam".to_string(),
            },
        ];
        let summary = build_summary(&genes);
        assert_eq!(summary.get("beta-lactam"), Some(&2));
        assert_eq!(summary.get("aminoglycoside"), Some(&1));
    }

    #[test]
    fn test_parse_amrfinder_tsv_mock() {
        let temp_dir = tempfile::TempDir::new().unwrap();
        let tsv_path = temp_dir.path().join("amr_output.tsv");
        let content = "gene_name\tscope\ttarget_type\tmethod\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tannotation\tproduct\tresistance\n\
TEM-1\tCORE\tGENE\tEXACT\t100.0\t100.0\tcontig1\t100\t200\t+\tTEM-1 beta-lactamase\tbeta-lactamase\tbeta-lactam\n\
aac(3)-IIa\tACCESSORY\tGENE\tALLELE\t99.5\t100.0\tcontig2\t300\t400\t+\taac(3)-IIa\taminoglycoside acetyltransferase\taminoglycoside\n";
        std::fs::write(&tsv_path, content).unwrap();
        let genes = parse_amrfinder_tsv(&tsv_path).unwrap();
        assert_eq!(genes.len(), 2);
        assert_eq!(genes[0].gene_name, "TEM-1");
        assert_eq!(genes[1].gene_name, "aac(3)-IIa");
    }

    #[test]
    fn test_write_and_summary() {
        let genes = vec![AmrGene {
            gene_name: "TEM-1".to_string(), scope: "CORE".to_string(),
            target_type: "GENE".to_string(), method: "EXACT".to_string(),
            identity: 100.0, coverage: 100.0, contig: "contig1".to_string(),
            start: 100, end: 200, strand: "+".to_string(),
            annotation: "TEM-1 beta-lactamase".to_string(),
            product: "beta-lactamase".to_string(), resistance: "beta-lactam".to_string(),
        }];
        let summary = build_summary(&genes);
        let result = AmrFinderResult { genes, summary_by_class: summary };
        let temp_dir = tempfile::TempDir::new().unwrap();
        result.write_to(temp_dir.path()).unwrap();
        assert!(temp_dir.path().join("amr_results.tsv").exists());
        assert!(temp_dir.path().join("amr_summary.txt").exists());
    }

    #[test]
    fn test_convert_csv_to_gff() {
        let temp_dir = tempfile::TempDir::new().unwrap();
        let csv_path = temp_dir.path().join("gene_data.csv");
        let gff_path = temp_dir.path().join("gene_data.gff");
        let csv_content = "cluster_id,annotation,product,contig,start,end,strand\ngene_001,hypothetical protein,hypothetical,contig1,100,200,+\ngene_002,transposase,transposase,contig1,300,400,-";
        std::fs::write(&csv_path, csv_content).unwrap();
        convert_csv_to_gff(&csv_path, &gff_path).unwrap();
        let gff_content = fs::read_to_string(&gff_path).unwrap();
        assert!(gff_content.contains("##gff-version 3"));
        assert!(gff_content.contains("gene_001"));
    }
}
