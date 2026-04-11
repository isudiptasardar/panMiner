//! Sourmash-based genome distance estimation.
//!
//! Computes pairwise ANI/distance matrices between genomes using
//! sourmash MinHash sketches. Feature-gated behind `sourmash` feature.
//!
//! When the `sourmash` feature is not enabled, falls back to logging
//! a warning and returning an empty result.

use std::path::PathBuf;

use crate::error::Result;

/// Result of an MDS projection.
pub struct MdsProjection {
    /// 2D coordinates for each genome.
    pub coordinates: Vec<(f64, f64)>,
    /// Genome labels.
    pub labels: Vec<String>,
}

/// Compute a pairwise distance matrix using sourmash.
///
/// When the `sourmash` feature is not enabled, this logs a warning
/// and returns an error.
pub fn compute_distance_matrix(_genomes: &[PathBuf]) -> Result<Vec<Vec<f64>>> {
    #[cfg(feature = "sourmash")]
    {
        compute_distance_matrix_sourmash(_genomes)
    }
    #[cfg(not(feature = "sourmash"))]
    {
        tracing::warn!("Sourmash feature not enabled. Install with --features sourmash for distance estimation.");
        Err(crate::Error::Output(
            "Sourmash not available. Enable the 'sourmash' feature or use FastANI instead.".to_string()
        ))
    }
}

/// Compute MDS projection from a distance matrix.
///
/// Uses classical MDS (metric MDS) to project genomes into 2D space
/// for visualization.
pub fn compute_mds(distances: &[Vec<f64>]) -> Result<MdsProjection> {
    let n = distances.len();
    if n == 0 {
        return Ok(MdsProjection {
            coordinates: Vec::new(),
            labels: Vec::new(),
        });
    }

    // Classical MDS: eigen-decomposition of centered distance matrix
    // This is a simplified implementation for 2D projection

    // Step 1: Square the distance matrix
    let mut d_sq = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            d_sq[i][j] = distances[i][j] * distances[i][j];
        }
    }

    // Step 2: Double centering: B = -0.5 * (I - 1/n * 11') * D^2 * (I - 1/n * 11')
    let mut b = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let row_mean: f64 = d_sq[i].iter().sum::<f64>() / n as f64;
            let col_mean: f64 = (0..n).map(|k| d_sq[k][j]).sum::<f64>() / n as f64;
            let grand_mean: f64 = d_sq.iter().flat_map(|r| r.iter()).sum::<f64>() / (n * n) as f64;
            b[i][j] = -0.5 * (d_sq[i][j] - row_mean - col_mean + grand_mean);
        }
    }

    // Step 3: Power iteration to find top 2 eigenvectors
    let mut coords = Vec::new();
    for dim in 0..2 {
        let mut v = vec![1.0 / (n as f64).sqrt(); n];
        for _ in 0..100 {
            // Multiply: v' = B * v
            let mut v_new = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    v_new[i] += b[i][j] * v[j];
                }
            }
            // Normalize
            let norm: f64 = v_new.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm > 0.0 {
                for x in v_new.iter_mut() {
                    *x /= norm;
                }
            }
            v = v_new;
        }
        // Eigenvalue approximation
        let mut bv = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                bv[i] += b[i][j] * v[j];
            }
        }
        let eigenvalue: f64 = v.iter().zip(bv.iter()).map(|(vi, bvi)| vi * bvi).sum();
        let sqrt_eigenvalue = if eigenvalue > 0.0 { eigenvalue.sqrt() } else { 0.0 };

        // Deflate: B = B - lambda * v * v'
        for i in 0..n {
            for j in 0..n {
                b[i][j] -= eigenvalue * v[i] * v[j];
            }
        }

        coords.push((v, sqrt_eigenvalue));
    }

    // Build 2D coordinates
    let mut coordinates = Vec::new();
    for i in 0..n {
        let x = coords[0].0[i] * coords[0].1;
        let y = coords[1].0[i] * coords[1].1;
        coordinates.push((x, y));
    }

    let labels = (0..n).map(|i| format!("genome_{}", i)).collect();

    Ok(MdsProjection {
        coordinates,
        labels,
    })
}

#[cfg(feature = "sourmash")]
fn compute_distance_matrix_sourmash(genomes: &[PathBuf]) -> Result<Vec<Vec<f64>>> {
    // When sourmash feature is enabled, use the sourmash crate
    // for accurate MinHash-based ANI computation.
    // For now, return a placeholder since the sourmash crate API
    // may vary between versions.
    tracing::info!("Computing distance matrix for {} genomes using sourmash", genomes.len());
    Err(crate::Error::Output(
        "Sourmash integration not yet implemented. Use FastANI instead.".to_string()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_distance_matrix_no_feature() {
        // When sourmash feature is not enabled, should return error
        let genomes = vec![PathBuf::from("test.fna")];
        let result = compute_distance_matrix(&genomes);
        #[cfg(not(feature = "sourmash"))]
        assert!(result.is_err());
    }

    #[test]
    fn test_compute_mds_empty() {
        let result = compute_mds(&[]).unwrap();
        assert!(result.coordinates.is_empty());
        assert!(result.labels.is_empty());
    }

    #[test]
    fn test_compute_mds_identity() {
        // Identity distance matrix (all zeros on diagonal)
        let distances = vec![
            vec![0.0, 0.05, 0.1],
            vec![0.05, 0.0, 0.08],
            vec![0.1, 0.08, 0.0],
        ];
        let result = compute_mds(&distances).unwrap();
        assert_eq!(result.coordinates.len(), 3);
        assert_eq!(result.labels.len(), 3);
    }

    #[test]
    fn test_mds_projection_labels() {
        let distances = vec![
            vec![0.0, 0.1],
            vec![0.1, 0.0],
        ];
        let result = compute_mds(&distances).unwrap();
        assert_eq!(result.labels.len(), 2);
        assert_eq!(result.labels[0], "genome_0");
        assert_eq!(result.labels[1], "genome_1");
    }
}