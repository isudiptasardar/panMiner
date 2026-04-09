//! GFF3 parser with memory-mapped file support.

use std::path::Path;
use std::collections::HashMap;
use rayon::prelude::*;

use crate::error::Result;
use crate::graph::{Gene, GeneId, GenomeId, Strand};
use super::mmap::MmapFile;

/// GFF3 parser with memory-mapped file support.
///
/// This parser reads GFF3 files efficiently using memory mapping
/// and can parse multiple files in parallel using Rayon.
pub struct GffParser {
    mmap: MmapFile,
    genome_id: GenomeId,
}

/// A parsed GFF record.
#[derive(Debug, Clone)]
pub struct GffRecord {
    pub seqid: String,
    pub feature_type: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    pub attributes: HashMap<String, String>,
}

impl GffParser {
    /// Open a GFF3 file for parsing.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the GFF3 file
    /// * `genome_id` - Identifier for this genome
    ///
    /// # Example
    ///
    /// ```no_run
    /// use panminer::io::GffParser;
    /// use panminer::graph::GenomeId;
    ///
    /// let parser = GffParser::open(std::path::Path::new("genome.gff"), GenomeId::new("sample1")).unwrap();
    /// let genes = parser.parse_genes().unwrap();
    /// ```
    pub fn open(path: &Path, genome_id: GenomeId) -> Result<Self> {
        let mmap = MmapFile::open(path)?;
        Ok(Self { mmap, genome_id })
    }

    /// Parse all gene features from the GFF3 file.
    ///
    /// This method extracts all features of type "gene" or "CDS"
    /// and returns them as Gene objects.
    pub fn parse_genes(&self) -> Result<Vec<Gene>> {
        let bytes = self.mmap.as_bytes();

        // Separate GFF lines and FASTA part
        let mut gff_bytes = bytes;
        let mut fasta_bytes: Option<&[u8]> = None;

        if let Some(pos) = bytes.windows(8).position(|w| w == b"##FASTA\n" || w == b"##FASTA\r") {
            gff_bytes = &bytes[..pos];
            // Find the start of the next line (the FASTA sequence)
            if let Some(next_line_pos) = bytes[pos..].iter().position(|&b| b == b'\n') {
                fasta_bytes = Some(&bytes[pos + next_line_pos + 1..]);
            }
        }

        // Find line boundaries for GFF part
        let lines: Vec<&[u8]> = gff_bytes.split(|&b| b == b'\n').collect();

        // Parse lines in parallel
        let mut genes: Vec<Gene> = lines
            .par_iter()
            .filter_map(|line| {
                // Skip comments and empty lines
                if line.starts_with(b"#") || line.is_empty() {
                    return None;
                }

                // Parse the GFF record
                let record = self.parse_line(line)?;
                let gene = record_to_gene(&record, &self.genome_id);
                gene
            })
            .collect();

        // If we have FASTA data, extract sequences
        if let Some(fasta) = fasta_bytes {
            let fasta_parser = crate::io::fasta::FastaIterator::new(fasta);
            let contigs: HashMap<String, crate::graph::Sequence> = fasta_parser
                .map(|record| (record.id, record.sequence))
                .collect();

            genes.par_iter_mut().for_each(|gene| {
                if let Some(contig_seq) = contigs.get(&gene.contig) {
                    if gene.start > 0 && gene.end <= contig_seq.len() && gene.start <= gene.end {
                        let mut seq = contig_seq[gene.start - 1..gene.end].to_vec();

                        // Reverse complement if negative strand
                        if gene.strand == Strand::Minus {
                            seq.reverse();
                            for byte in seq.iter_mut() {
                                *byte = match *byte {
                                    b'A' | b'a' => b'T',
                                    b'C' | b'c' => b'G',
                                    b'G' | b'g' => b'C',
                                    b'T' | b't' => b'A',
                                    _ => *byte,
                                };
                            }
                        }

                        gene.sequence = seq;
                    }
                }
            });
        }

        Ok(genes)
    }

    /// Parse a single GFF3 line into a record.
    fn parse_line(&self, line: &[u8]) -> Option<GffRecord> {
        // Strip trailing \r (Windows line endings)
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        let fields: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();

        if fields.len() < 9 {
            return None;
        }

        let seqid = String::from_utf8_lossy(fields[0]).to_string();
        let feature_type = String::from_utf8_lossy(fields[2]).to_string();

        let start = String::from_utf8_lossy(fields[3])
            .parse::<usize>()
            .ok()?;
        let end = String::from_utf8_lossy(fields[4])
            .parse::<usize>()
            .ok()?;

        let strand = parse_strand(fields[6])?;

        let attributes = parse_attributes(fields[8]);

        Some(GffRecord {
            seqid,
            feature_type,
            start,
            end,
            strand,
            attributes,
        })
    }

    /// Get the genome ID for this file.
    pub fn genome_id(&self) -> &GenomeId {
        &self.genome_id
    }
}

/// Parse the strand field.
fn parse_strand(s: &[u8]) -> Option<Strand> {
    match s {
        b"+" => Some(Strand::Plus),
        b"-" => Some(Strand::Minus),
        b"." => Some(Strand::Unknown),
        _ => None,
    }
}

/// Parse GFF3 attributes (column 9).
fn parse_attributes(attrs: &[u8]) -> HashMap<String, String> {
    let mut map = HashMap::new();

    for pair in attrs.split(|&b| b == b';') {
        if pair.is_empty() {
            continue;
        }

        if let Some(eq_pos) = pair.iter().position(|&b| b == b'=') {
            let key = String::from_utf8_lossy(&pair[..eq_pos]).to_string();
            let value = String::from_utf8_lossy(&pair[eq_pos + 1..]).to_string();
            map.insert(key, value);
        }
    }

    map
}

/// Convert a GFF record to a Gene.
fn record_to_gene(record: &GffRecord, genome_id: &GenomeId) -> Option<Gene> {
    // Only process gene and CDS features
    if record.feature_type != "gene" && record.feature_type != "CDS" {
        return None;
    }

    // Get the gene ID from attributes
    let gene_id = record.attributes
        .get("ID")
        .or_else(|| record.attributes.get("Name"))
        .cloned()
        .unwrap_or_else(|| format!("{}_{}_{}", record.seqid, record.start, record.end));

    // Get annotation
    let annotation = record.attributes
        .get("product")
        .or_else(|| record.attributes.get("gene"))
        .cloned();

    Some(Gene {
        id: GeneId::new(gene_id),
        sequence: Vec::new(), // Sequence loaded separately from FASTA
        genome_id: genome_id.clone(),
        contig: record.seqid.clone(),
        start: record.start,
        end: record.end,
        strand: record.strand,
        annotation,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_gff() -> Result<NamedTempFile> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "##gff-version 3")?;
        writeln!(temp, "seq1\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene1;product=test gene")?;
        writeln!(temp, "seq1\tProkka\tCDS\t100\t200\t.\t+\t0\tID=cds1;Parent=gene1")?;
        writeln!(temp, "seq2\tProkka\tgene\t500\t800\t.\t-\t.\tID=gene2;product=another gene")?;
        Ok(temp)
    }

    fn create_test_gff_with_fasta() -> Result<NamedTempFile> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "##gff-version 3")?;
        writeln!(temp, "seq1\tProkka\tgene\t1\t4\t.\t+\t.\tID=gene1")?;
        writeln!(temp, "seq2\tProkka\tgene\t2\t5\t.\t-\t.\tID=gene2")?;
        writeln!(temp, "##FASTA")?;
        writeln!(temp, ">seq1")?;
        writeln!(temp, "ATCGATCG")?;
        writeln!(temp, ">seq2")?;
        writeln!(temp, "GCTAGCTA")?;
        Ok(temp)
    }

    #[test]
    fn test_parse_strand() {
        assert_eq!(parse_strand(b"+"), Some(Strand::Plus));
        assert_eq!(parse_strand(b"-"), Some(Strand::Minus));
        assert_eq!(parse_strand(b"."), Some(Strand::Unknown));
        assert_eq!(parse_strand(b"x"), None);
    }

    #[test]
    fn test_parse_attributes() {
        let attrs = b"ID=gene1;product=test gene;Name=gene_1";
        let map = parse_attributes(attrs);
        assert_eq!(map.get("ID"), Some(&"gene1".to_string()));
        assert_eq!(map.get("product"), Some(&"test gene".to_string()));
        assert_eq!(map.get("Name"), Some(&"gene_1".to_string()));
    }

    #[test]
    fn test_gff_parser() -> Result<()> {
        let temp = create_test_gff()?;
        let parser = GffParser::open(temp.path(), GenomeId::new("test_genome"))?;
        let genes = parser.parse_genes().unwrap();

        // Should have 3 gene/CDS features
        assert_eq!(genes.len(), 3);

        Ok(())
    }

    #[test]
    fn test_gff_parser_with_fasta() -> Result<()> {
        let temp = create_test_gff_with_fasta()?;
        let parser = GffParser::open(temp.path(), GenomeId::new("test_genome"))?;
        let genes = parser.parse_genes().unwrap();

        assert_eq!(genes.len(), 2);

        // Sequence should be extracted
        assert_eq!(genes[0].sequence, b"ATCG".to_vec()); // seq1: 1-4

        // Negative strand should be reverse complemented in a real implementation,
        // but for now let's just ensure it's extracted (CTAG)
        assert_eq!(genes[1].sequence, b"CTAG".to_vec()); // seq2: 2-5

        Ok(())
    }
}