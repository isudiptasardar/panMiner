//! DNA-to-protein translation using the standard genetic code (Table 11).
//!
//! Provides translation for bacterial gene sequences, including
//! standard codons, start codon overrides, and stop codon handling.

use std::sync::OnceLock;

/// Start codons that translate to Methionine (M) instead of their
/// standard amino acid. Table 11 bacterial start codons.
const START_CODONS: &[[u8; 3]] = &[*b"ATG", *b"GTG", *b"TTG"];

/// Cached codon table (built once on first use).
static CODON_TABLE: OnceLock<[u8; 64]> = OnceLock::new();

/// Get the codon translation table (cached after first call).
fn codon_table() -> &'static [u8; 64] {
    CODON_TABLE.get_or_init(build_codon_table)
}

/// Translate a DNA sequence to amino acids using the standard genetic code.
///
/// - Stops at the first stop codon (`*`)
/// - Handles partial codons by truncating (does not pad)
/// - Handles start codons: ATG/GTG/TTG → M at position 0
/// - Handles ambiguous bases (N, R, Y, etc.) → X
pub fn translate(dna: &[u8]) -> Vec<u8> {
    if dna.len() < 3 {
        return Vec::new();
    }

    let mut protein = Vec::with_capacity(dna.len() / 3);
    let table = codon_table();

    for i in (0..dna.len() - 2).step_by(3) {
        let codon = &dna[i..i + 3];
        let aa = if i == 0 && is_start_codon(codon) {
            b'M'
        } else {
            lookup_codon(codon, table)
        };

        if aa == b'*' {
            break;
        }
        protein.push(aa);
    }

    protein
}

/// Translate DNA to protein, returning the full sequence including
/// the stop codon marker. Used for checking for internal stop codons.
pub fn translate_with_stop(dna: &[u8]) -> Vec<u8> {
    if dna.len() < 3 {
        return Vec::new();
    }

    let mut protein = Vec::with_capacity(dna.len() / 3);
    let table = codon_table();

    for i in (0..dna.len() - 2).step_by(3) {
        let codon = &dna[i..i + 3];
        let aa = if i == 0 && is_start_codon(codon) {
            b'M'
        } else {
            lookup_codon(codon, table)
        };
        protein.push(aa);
    }

    protein
}

/// Check if a codon is a start codon using direct byte comparison (no allocation).
fn is_start_codon(codon: &[u8]) -> bool {
    if codon.len() != 3 {
        return false;
    }
    START_CODONS.iter().any(|sc| {
        codon[0].to_ascii_uppercase() == sc[0]
            && codon[1].to_ascii_uppercase() == sc[1]
            && codon[2].to_ascii_uppercase() == sc[2]
    })
}

/// Build a lookup table mapping codon indices to amino acids.
fn build_codon_table() -> [u8; 64] {
    let mut table = [b'X'; 64];
    // T = 0, C = 1, A = 2, G = 3
    // Index = first*16 + second*4 + third
    let codons = [
        "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
        "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
        "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
        "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
        "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
        "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
        "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
        "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG",
    ];
    let aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

    for (i, codon) in codons.iter().enumerate() {
        if let Some(idx) = codon_to_index(codon.as_bytes()) {
            table[idx] = aas.as_bytes()[i];
        }
    }
    table
}

/// Convert a codon to a lookup table index.
fn codon_to_index(codon: &[u8]) -> Option<usize> {
    if codon.len() != 3 {
        return None;
    }
    let n1 = nuc_to_index(codon[0])?;
    let n2 = nuc_to_index(codon[1])?;
    let n3 = nuc_to_index(codon[2])?;
    Some(n1 * 16 + n2 * 4 + n3)
}

/// Convert a nucleotide to its index (T=0, C=1, A=2, G=3).
fn nuc_to_index(n: u8) -> Option<usize> {
    match n.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Look up a codon in the translation table.
fn lookup_codon(codon: &[u8], table: &[u8; 64]) -> u8 {
    match codon_to_index(codon) {
        Some(idx) => table[idx],
        None => b'X',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_basic() {
        let dna = b"ATGTTTTTCTAA";
        let protein = translate(dna);
        assert_eq!(protein, b"MFF");
    }

    #[test]
    fn test_translate_start_codon_gtg() {
        let dna = b"GTGTTTTTCTAA";
        let protein = translate(dna);
        assert_eq!(protein[0], b'M');
    }

    #[test]
    fn test_translate_start_codon_ttg() {
        let dna = b"TTGTTTTTCTAA";
        let protein = translate(dna);
        assert_eq!(protein[0], b'M');
    }

    #[test]
    fn test_translate_stop_early() {
        let dna = b"ATGTAA";
        let protein = translate(dna);
        assert_eq!(protein, b"M");
    }

    #[test]
    fn test_translate_with_stop() {
        let dna = b"ATGTTTTTCTAA";
        let protein = translate_with_stop(dna);
        assert_eq!(protein, b"MFF*");
    }

    #[test]
    fn test_translate_short_sequence() {
        let dna = b"AT";
        let protein = translate(dna);
        assert!(protein.is_empty());
    }

    #[test]
    fn test_translate_ambiguous_base() {
        let dna = b"ATGNTTTTCTAA";
        let protein = translate(dna);
        assert_eq!(protein[0], b'M');
        assert_eq!(protein[1], b'X');
    }

    #[test]
    fn test_translate_all_stops() {
        let dna1 = b"ATGTTTTTCTAA";
        assert_eq!(translate(dna1), b"MFF");

        let dna2 = b"ATGTTTTTCTAG";
        assert_eq!(translate(dna2), b"MFF");

        let dna3 = b"ATGTTTTTCTGA";
        assert_eq!(translate(dna3), b"MFF");
    }

    #[test]
    fn test_translate_lowercase() {
        let dna = b"atgtttttctaa";
        let protein = translate(dna);
        assert_eq!(protein, b"MFF");
    }

    #[test]
    fn test_codon_table_cached() {
        let t1 = codon_table();
        let t2 = codon_table();
        assert!(std::ptr::eq(t1, t2), "codon_table should be cached");
    }
}