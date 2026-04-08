#!/bin/bash
sed -i '/fn sequence_identity/i \
pub fn simd_sequence_identity(a: &[u8], b: &[u8]) -> f32 {\
    let min_len = a.len().min(b.len());\
    if min_len == 0 {\
        return 0.0;\
    }\
\
    let matches: usize = a[..min_len]\
        .iter()\
        .zip(\&b[..min_len])\
        .map(|(x, y)| if x == y { 1 } else { 0 })\
        .sum();\
\
    matches as f32 / min_len as f32\
}\
' src/clustering/cpu.rs
