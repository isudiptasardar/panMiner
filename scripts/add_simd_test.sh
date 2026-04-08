#!/bin/bash
cat << 'INNER_EOF' >> src/clustering/cpu.rs

#[cfg(test)]
mod simd_tests {
    use super::*;

    #[test]
    fn test_simd_sequence_identity() {
        let seq_a = b"ATCGATCGATCGATCGATCGATCGATCGATCG";
        let seq_b = b"ATCGATCGATCGATCGATCGATCGATCGATCG";
        let seq_c = b"ATCGATCGATCGATCGATCGATCGATCGAACG"; // 1 mismatch
        
        assert_eq!(simd_sequence_identity(seq_a, seq_b), 1.0);
        assert_eq!(simd_sequence_identity(seq_a, seq_c), 31.0 / 32.0);
        
        // Empty
        assert_eq!(simd_sequence_identity(b"", b"A"), 0.0);
        
        // Different lengths
        assert_eq!(simd_sequence_identity(b"ATCG", b"ATC"), 1.0); // limited by min_len
    }
}
INNER_EOF
