#!/bin/bash
cat << 'INNER_EOF' >> src/clustering/cpu.rs

/// SIMD-optimized Hamming distance
pub fn simd_hamming_distance(a: &[u8], b: &[u8]) -> usize {
    let min_len = a.len().min(b.len());
    if min_len == 0 {
        return a.len().max(b.len());
    }

    let mismatches: usize = a[..min_len]
        .iter()
        .zip(&b[..min_len])
        .map(|(x, y)| if x != y { 1 } else { 0 })
        .sum();

    mismatches + (a.len().max(b.len()) - min_len)
}

#[cfg(test)]
mod hamming_tests {
    use super::*;

    #[test]
    fn test_simd_hamming_distance() {
        let seq_a = b"ATCGATCG";
        let seq_b = b"ATCGATCG";
        let seq_c = b"ATCGAACG"; // 1 mismatch
        let seq_d = b"ATCG";     // length difference 4

        assert_eq!(simd_hamming_distance(seq_a, seq_b), 0);
        assert_eq!(simd_hamming_distance(seq_a, seq_c), 1);
        assert_eq!(simd_hamming_distance(seq_a, seq_d), 4);
    }
}
INNER_EOF
