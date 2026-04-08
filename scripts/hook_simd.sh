#!/bin/bash
sed -i 's/self.search_sequence(query_seq, contig_sequences)/self.simd_search_sequence(query_seq, contig_sequences)/g' src/correction/missing.rs
sed -i 's/Self::sequence_identity/Self::simd_sequence_identity/g' src/clustering/cpu.rs
