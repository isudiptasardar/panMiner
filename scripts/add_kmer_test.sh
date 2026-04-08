#!/bin/bash
cat << 'INNER_EOF' >> src/correction/missing.rs

#[cfg(test)]
mod simd_tests {
    use super::*;

    #[test]
    fn test_simd_kmer_search() {
        let recoverer = MissingGeneRecoverer::new();
        let query = b"ATCGATCGATCGATCG";
        let mut contigs = HashMap::new();
        contigs.insert("contig1".to_string(), b"NNNATCGATCGATCGATCGNNN".to_vec());
        
        assert!(recoverer.simd_search_sequence(query, &contigs));
        
        let mut contigs_fail = HashMap::new();
        contigs_fail.insert("contig1".to_string(), b"NNNATCGAACGATCGAACGNNN".to_vec());
        assert!(!recoverer.simd_search_sequence(query, &contigs_fail));
    }
}
INNER_EOF
