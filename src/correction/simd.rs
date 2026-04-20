//! SIMD sequence comparison utilities.
//!
//! Uses runtime CPU feature detection for accelerated sequence comparison.
//! Falls back to scalar implementation on unsupported architectures.
//! Also provides edlib-based Levenshtein alignment for Panaroo compatibility.

#[cfg(target_arch = "x86_64")]
use std::arch::is_x86_feature_detected;

/// Calculate sequence identity using the best available method.
///
/// Uses SIMD (AVX2/NEON) when available, falls back to scalar comparison.
/// Returns the fraction of matching positions between two sequences.
pub fn compare_sequences(a: &[u8], b: &[u8]) -> f64 {
    let len_a = a.len();
    let len_b = b.len();
    let min_len = len_a.min(len_b);

    if min_len == 0 {
        return 0.0;
    }

    // Try AVX2 first (x86_64)
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx2") {
        return simd_sequence_identity_avx2(a, b);
    }

    // Try NEON (aarch64)
    #[cfg(target_arch = "aarch64")]
    if is_aarch64_feature_detected!("neon") {
        return simd_sequence_identity_neon(a, b);
    }

    // Fallback to scalar comparison
    scalar_sequence_identity(a, b)
}

/// Scalar sequence identity calculation (fallback).
fn scalar_sequence_identity(a: &[u8], b: &[u8]) -> f64 {
    let min_len = a.len().min(b.len());
    let matches: usize = a.iter().zip(b).filter(|(a, b)| a == b).count();
    matches as f64 / min_len as f64
}

/// AVX2-accelerated sequence comparison.
#[cfg(target_arch = "x86_64")]
fn simd_sequence_identity_avx2(a: &[u8], b: &[u8]) -> f64 {
    use std::arch::x86_64::*;

    let min_len = a.len().min(b.len());
    let mut matches = 0u64;
    let mut i = 0;

    // Process 32 bytes at a time using AVX2
    while i + 31 < min_len {
        unsafe {
            let va = _mm256_loadu_si256(a.as_ptr().add(i) as *const _);
            let vb = _mm256_loadu_si256(b.as_ptr().add(i) as *const _);
            let eq = _mm256_cmpeq_epi8(va, vb);
            let mask = _mm256_movemask_epi8(eq);
            matches += mask.count_ones() as u64;
        }
        i += 32;
    }

    // Handle remaining bytes with scalar
    while i < min_len {
        if a[i] == b[i] {
            matches += 1;
        }
        i += 1;
    }

    matches as f64 / min_len as f64
}

/// NEON-accelerated sequence comparison (ARM).
#[cfg(target_arch = "aarch64")]
fn simd_sequence_identity_neon(a: &[u8], b: &[u8]) -> f64 {
    use std::arch::aarch64::*;

    let min_len = a.len().min(b.len());
    let mut matches = 0u64;
    let mut i = 0;

    // Process 16 bytes at a time using NEON
    while i + 15 < min_len {
        unsafe {
            let va = vld1q_u8(a.as_ptr().add(i));
            let vb = vld1q_u8(b.as_ptr().add(i));
            let eq = vceqq_u8(va, vb);
            // vceqq_u8 returns 0xFF for matching bytes, 0x00 otherwise.
            // Use horizontal add to count all 16 matching bytes.
            matches += vaddvq_u8(eq) as u64 / 255;
        }
        i += 16;
    }

    // Handle remaining bytes with scalar
    while i < min_len {
        if a[i] == b[i] {
            matches += 1;
        }
        i += 1;
    }

    matches as f64 / min_len as f64
}

/// Compute Levenshtein (edit) distance between two sequences.
/// Returns the normalized identity (1.0 - edit_distance / alignment_length).
///
/// Uses a standard dynamic programming approach for global alignment.
/// For semi-global alignment (query within target), use `align_semiglobal`.
pub fn align_sequences(query: &[u8], target: &[u8]) -> f64 {
    if query.is_empty() || target.is_empty() {
        return 0.0;
    }

    let distance = levenshtein_distance(query, target);
    let alignment_length = query.len().max(target.len());

    1.0 - (distance as f64 / alignment_length as f64)
}

/// Compute semi-global alignment of a query within a target.
/// Returns (identity, edit_distance, alignment_length).
///
/// Uses HW (semi-global/infix) mode: the query is aligned within the target,
/// allowing free gaps at the start and end of the target. This matches
/// Panaroo's find_missing approach using edlib with mode="HW".
pub fn align_semiglobal(query: &[u8], target: &[u8]) -> (f64, i32, usize) {
    if query.is_empty() || target.is_empty() {
        return (0.0, -1, 0);
    }

    let distance = levenshtein_distance_hw(query, target);
    let align_len = query.len();
    let identity = 1.0 - (distance as f64 / align_len as f64);

    (identity, distance, align_len)
}

/// Compute Levenshtein edit distance (global alignment).
///
/// Standard O(n*m) dynamic programming implementation.
fn levenshtein_distance(a: &[u8], b: &[u8]) -> i32 {
    let len_a = a.len();
    let len_b = b.len();

    if len_a == 0 { return len_b as i32; }
    if len_b == 0 { return len_a as i32; }

    // Use only two rows for space efficiency
    let mut prev_row: Vec<i32> = (0..=len_b as i32).collect();
    let mut curr_row: Vec<i32> = vec![0; len_b + 1];

    for (i, &a_byte) in a.iter().enumerate() {
        curr_row[0] = (i + 1) as i32;
        for (j, &b_byte) in b.iter().enumerate() {
            let cost = if a_byte == b_byte { 0 } else { 1 };
            curr_row[j + 1] = (prev_row[j + 1] + 1)      // deletion
                .min(curr_row[j] + 1)                      // insertion
                .min(prev_row[j] + cost);                  // substitution
        }
        std::mem::swap(&mut prev_row, &mut curr_row);
    }

    prev_row[len_b]
}

/// Compute Levenshtein edit distance in HW (semi-global/infix) mode.
///
/// In HW mode, the query is aligned within the target, allowing free
/// leading and trailing gaps in the target. This finds the best match
/// of the query anywhere within the target sequence.
///
/// This matches edlib's HW mode used by Panaroo for find_missing.
fn levenshtein_distance_hw(query: &[u8], target: &[u8]) -> i32 {
    let query_len = query.len();
    let target_len = target.len();

    if query_len == 0 { return 0; }
    if target_len == 0 { return query_len as i32; }

    // Use two rows for space efficiency
    let mut prev_row: Vec<i32> = (0..=query_len as i32).collect();
    let mut curr_row: Vec<i32> = vec![0; query_len + 1];

    // HW mode: free leading gaps in target (start curr_row[0] = 0 is implicit)
    // Free trailing gaps: take the minimum of the last row
    let mut best_distance = query_len as i32;

    for (_j, &t_byte) in target.iter().enumerate() {
        curr_row[0] = 0; // Free leading gaps in target (HW mode)
        for (i, &q_byte) in query.iter().enumerate() {
            let cost = if q_byte == t_byte { 0 } else { 1 };
            curr_row[i + 1] = (prev_row[i + 1] + 1)       // deletion from query
                .min(curr_row[i] + 1)                       // insertion into query
                .min(prev_row[i] + cost);                   // substitution/match
        }
        // Track the best (minimum) distance at the end of the target row
        best_distance = best_distance.min(curr_row[query_len]);
        std::mem::swap(&mut prev_row, &mut curr_row);
    }

    best_distance
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_sequences() {
        let a = b"ATCGATCGATCG";
        let b = b"ATCGATCGATCG";
        let identity = compare_sequences(a, b);
        assert!((identity - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_different_sequences() {
        // ATCGATCGATCG vs GGGGGGGGGGGG: only positions 3, 7, 11 have G==G
        // 3/12 = 0.25
        let a = b"ATCGATCGATCG";
        let b = b"GGGGGGGGGGGG";
        let identity = compare_sequences(a, b);
        assert_eq!(identity, 0.25, "Expected 3/12 = 0.25 identity");
    }

    #[test]
    fn test_partial_match() {
        // ATCGATCG vs ATCGGGGG: positions 0,1,2,3,7 match (A,A,A,A,G)
        // 5/8 = 0.625
        let a = b"ATCGATCG";
        let b = b"ATCGGGGG";
        let identity = compare_sequences(a, b);
        assert_eq!(identity, 0.625, "Expected 5/8 = 0.625 identity");
    }

    #[test]
    fn test_empty_sequences() {
        let a = b"";
        let b = b"ATCG";
        let identity = compare_sequences(a, b);
        assert_eq!(identity, 0.0);
    }

    #[test]
    fn test_different_lengths() {
        let a = b"ATCGATCG";
        let b = b"ATCG";
        let identity = compare_sequences(a, b);
        // Should compare up to shorter length (4 matches out of 4)
        assert!((identity - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_longer_sequence() {
        let a = b"ATCG".repeat(100);
        let b = b"ATCG".repeat(100);
        let identity = compare_sequences(&a, &b);
        assert!((identity - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_levenshtein_distance_identical() {
        let dist = levenshtein_distance(b"ATCGATCG", b"ATCGATCG");
        assert_eq!(dist, 0);
    }

    #[test]
    fn test_levenshtein_distance_one_sub() {
        let dist = levenshtein_distance(b"ATCGATCG", b"ATCGATCC");
        assert_eq!(dist, 1);
    }

    #[test]
    fn test_align_sequences_identical() {
        let identity = align_sequences(b"ATCGATCG", b"ATCGATCG");
        assert!((identity - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_align_sequences_different() {
        let identity = align_sequences(b"ATCGATCG", b"GCTAGCTA");
        // 6 substitutions out of 8 = 0.25 identity
        assert!(identity < 0.5, "Expected low identity for completely different sequences, got {}", identity);
    }

    #[test]
    fn test_align_semiglobal_exact_match() {
        // Query found exactly within target
        let (identity, dist, len) = align_semiglobal(b"ATCG", b"XXXATCGYYY");
        assert_eq!(dist, 0);
        assert_eq!(len, 4);
        assert!((identity - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_align_semiglobal_near_match() {
        // Query found with one substitution within target
        let (identity, dist, _len) = align_semiglobal(b"ATCG", b"XXXATCAYYY");
        assert_eq!(dist, 1); // 1 substitution
        assert!((identity - 0.75).abs() < 0.001);
    }

    #[test]
    fn test_levenshtein_hw_embedded() {
        // Query "ABC" embedded in larger target "XYZABCPQR"
        // HW mode should find it with distance 0
        let dist = levenshtein_distance_hw(b"ABC", b"XYZABCPQR");
        assert_eq!(dist, 0, "HW mode should find exact embedded match");
    }

    #[test]
    fn test_levenshtein_hw_no_match() {
        // Query has no close match in target
        let dist = levenshtein_distance_hw(b"AAA", b"TTTTTT");
        assert_eq!(dist, 3, "No match should give edit distance = query length");
    }
}
