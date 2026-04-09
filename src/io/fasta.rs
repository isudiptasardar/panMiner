//! FASTA parser with memory-mapped file support.

use std::collections::HashMap;
use std::path::Path;

use crate::error::Result;
use crate::graph::Sequence;
use super::mmap::MmapFile;

/// FASTA parser with memory-mapped file support.
///
/// This parser reads FASTA files efficiently using memory mapping
/// and provides lazy iteration over sequences.
pub struct FastaParser {
    mmap: MmapFile,
}

/// A FASTA record containing a header and sequence.
#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// Header line (without the leading '>')
    pub header: String,
    /// Sequence identifier (first word of header)
    pub id: String,
    /// Sequence data
    pub sequence: Sequence,
}

impl FastaParser {
    /// Open a FASTA file for parsing.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the FASTA file
    ///
    /// # Example
    ///
    /// ```no_run
    /// use panminer::io::FastaParser;
    ///
    /// let parser = FastaParser::open(std::path::Path::new("sequences.fasta")).unwrap();
    /// let sequences = parser.parse_all().unwrap();
    /// ```
    pub fn open(path: &Path) -> Result<Self> {
        let mmap = MmapFile::open(path)?;
        Ok(Self { mmap })
    }

    /// Parse all sequences from the FASTA file.
    ///
    /// Returns a HashMap mapping sequence IDs to sequences.
    pub fn parse_all(&self) -> Result<HashMap<String, Sequence>> {
        let bytes = self.mmap.as_bytes();
        let mut sequences = HashMap::new();

        let mut current_header: Option<String> = None;
        let mut current_id: Option<String> = None;
        let mut current_seq: Sequence = Vec::new();

        for line in bytes.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }

            if line.starts_with(b">") {
                // Save previous record if exists
                if let (Some(_header), Some(id)) = (current_header.take(), current_id.take()) {
                    if !current_seq.is_empty() {
                        sequences.insert(id, current_seq.clone());
                        current_seq.clear();
                    }
                }

                // Parse header
                let header_line = &line[1..]; // Skip '>'
                let header_str = String::from_utf8_lossy(header_line).to_string();
                let id = header_str.split_whitespace().next().unwrap_or("").to_string();

                current_header = Some(header_str);
                current_id = Some(id);
            } else {
                // Append sequence line
                let line_str = String::from_utf8_lossy(line);
                let seq_bytes: Vec<u8> = line_str
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .map(|c| c as u8)
                    .collect();
                current_seq.extend(seq_bytes);
            }
        }

        // Don't forget the last record
        if let Some(id) = current_id {
            if !current_seq.is_empty() {
                sequences.insert(id, current_seq);
            }
        }

        Ok(sequences)
    }

    /// Iterate over sequences lazily.
    ///
    /// This is more memory-efficient for large files as it doesn't
    /// load all sequences into memory at once.
    pub fn iter(&self) -> FastaIterator<'_> {
        FastaIterator::new(self.mmap.as_bytes())
    }

    /// Get the file size in bytes.
    pub fn file_size(&self) -> usize {
        self.mmap.len()
    }
}

/// Lazy iterator over FASTA records.
pub struct FastaIterator<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> FastaIterator<'a> {
    pub fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, position: 0 }
    }
}

impl<'a> Iterator for FastaIterator<'a> {
    type Item = FastaRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position >= self.bytes.len() {
            return None;
        }

        // Find next header
        let remaining = &self.bytes[self.position..];

        // Skip to next '>'
        let header_start = remaining.iter().position(|&b| b == b'>')?;
        let header_end = header_start + remaining[header_start..]
            .iter()
            .position(|&b| b == b'\n')
            .unwrap_or(remaining.len() - header_start);

        let header_line = &remaining[header_start + 1..header_end];
        let header_str = String::from_utf8_lossy(header_line).to_string();
        let id = header_str.split_whitespace().next().unwrap_or("").to_string();

        // Find next header or end
        let seq_start = header_end + 1;
        let seq_end = remaining[seq_start..]
            .iter()
            .position(|&b| b == b'>')
            .map(|pos| seq_start + pos)
            .unwrap_or(remaining.len());

        // Extract sequence
        let seq_bytes: Sequence = remaining[seq_start..seq_end]
            .iter()
            .filter(|&&b| !b.is_ascii_whitespace())
            .copied()
            .collect();

        self.position = self.position + seq_end;

        Some(FastaRecord {
            header: header_str,
            id,
            sequence: seq_bytes,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_fasta() -> Result<NamedTempFile> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, ">seq1 description here")?;
        writeln!(temp, "ATCGATCG")?;
        writeln!(temp, ">seq2 another sequence")?;
        writeln!(temp, "GCTAGCTA")?;
        writeln!(temp, "NNNN")?;
        Ok(temp)
    }

    #[test]
    fn test_fasta_parser() -> Result<()> {
        let temp = create_test_fasta()?;
        let parser = FastaParser::open(temp.path())?;
        let sequences = parser.parse_all().unwrap();

        assert_eq!(sequences.len(), 2);
        assert_eq!(sequences.get("seq1"), Some(&b"ATCGATCG".to_vec()));
        assert_eq!(sequences.get("seq2"), Some(&b"GCTAGCTANNNN".to_vec()));

        Ok(())
    }

    #[test]
    fn test_fasta_iter() -> Result<()> {
        let temp = create_test_fasta()?;
        let parser = FastaParser::open(temp.path())?;

        let records: Vec<FastaRecord> = parser.iter().collect();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[1].id, "seq2");

        Ok(())
    }
}