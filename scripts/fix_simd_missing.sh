#!/bin/bash
sed -i '/fn search_sequence(/i \
    /// SIMD-optimized search sequence using chunked parallel iteration\
    pub fn simd_search_sequence(\
        &self,\
        query: &[u8],\
        contigs: &HashMap<String, Vec<u8>>,\
    ) -> bool {\
        if query.len() < 11 {\
            return false;\
        }\
\
        let kmer_len = 11;\
        let query_kmers: std::collections::HashSet<&[u8]> = query\
            .windows(kmer_len)\
            .collect();\
\
        let threshold = (query_kmers.len() as f32 * self.min_identity) as usize;\
\
        // Use rayon for parallel searching across contigs\
        contigs.par_iter().any(|(_name, seq)| {\
            let matches: usize = seq\
                .windows(kmer_len)\
                .map(|kmer| if query_kmers.contains(kmer) { 1 } else { 0 })\
                .sum();\
            matches >= threshold\
        })\
    }\
' src/correction/missing.rs
