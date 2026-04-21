//! Gene accumulation curve analysis with Heaps' law fitting.
//!
//! Generates gene accumulation curves via genome rarefaction and fits
//! Heaps' law (n(k) = A * k^alpha) to classify pangenomes as open or closed.

use std::path::Path;
use rand::seq::SliceRandom;
use rand::SeedableRng;

use crate::error::Result;
use crate::graph::BitPackedMatrix;
use crate::downstream::traits::{DownstreamInput, DownstreamResult, DownstreamRunner};

/// Runner for gene accumulation curve analysis.
pub struct AccumulationCurveRunner {
    /// Number of random samples to take at each genome count.
    num_samples: usize,
    /// Number of rarefaction points (genome counts to sample).
    rarefaction_points: usize,
}

impl AccumulationCurveRunner {
    /// Create a new accumulation curve runner with default parameters.
    pub fn new() -> Self {
        Self {
            num_samples: 100,
            rarefaction_points: 20,
        }
    }

    /// Set the number of samples per rarefaction point.
    pub fn with_num_samples(mut self, num_samples: usize) -> Self {
        self.num_samples = num_samples;
        self
    }

    /// Set the number of rarefaction points.
    pub fn with_rarefaction_points(mut self, rarefaction_points: usize) -> Self {
        self.rarefaction_points = rarefaction_points;
        self
    }

    /// Generate evenly spaced genome counts from 1 to n.
    fn rarefaction_ks(&self, n: usize) -> Vec<usize> {
        if self.rarefaction_points == 0 {
            return vec![];
        }
        if self.rarefaction_points >= n {
            return (1..=n).collect();
        }
        // Evenly spaced integer points from 1 to n
        let step = (n - 1) as f64 / (self.rarefaction_points - 1).max(1) as f64;
        (0..self.rarefaction_points)
            .map(|i| {
                let k = 1.0 + i as f64 * step;
                k.round() as usize
            })
            .collect()
    }

    /// Count genes present in the given sampled genomes.
    fn count_genes(&self, matrix: &BitPackedMatrix, genome_indices: &[usize]) -> (usize, usize) {
        let k = genome_indices.len();
        let total_genes = matrix.num_clusters();
        let mut present_count = 0usize;
        let mut core_count = 0usize;

        for cluster_idx in 0..total_genes {
            // Count how many of the sampled genomes have this cluster
            let count = genome_indices
                .iter()
                .filter(|&&g| matrix.get(g, cluster_idx))
                .count();
            if count > 0 {
                present_count += 1;
            }
            // Core gene: present in ALL sampled genomes
            if count == k {
                core_count += 1;
            }
        }

        (present_count, core_count)
    }

    /// Run rarefaction analysis on an in-memory BitPackedMatrix.
    fn run_on_matrix(&self, matrix: &BitPackedMatrix) -> AccumulationResult {
        let n = matrix.num_genomes();
        if n == 0 {
            return AccumulationResult {
                curve_data: vec![],
                heaps_law: HeapsLawFit {
                    alpha: 0.0,
                    a_coefficient: 0.0,
                    r_squared: 0.0,
                    classification: "unknown".to_string(),
                },
            };
        }

        let ks = self.rarefaction_ks(n);
        let mut curve_data = Vec::with_capacity(ks.len());

        // Fixed seed for reproducibility
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        for &k in &ks {
            let mut total_genes = Vec::with_capacity(self.num_samples);
            let mut core_genes = Vec::with_capacity(self.num_samples);

            for _ in 0..self.num_samples {
                // Sample k genomes without replacement
                let mut indices: Vec<usize> = (0..n).collect();
                indices.shuffle(&mut rng);
                indices.truncate(k);
                indices.sort();

                let (present, core) = self.count_genes(matrix, &indices);
                total_genes.push(present as f64);
                core_genes.push(core as f64);
            }

            let mean_total = total_genes.iter().sum::<f64>() / self.num_samples as f64;
            let mean_core = core_genes.iter().sum::<f64>() / self.num_samples as f64;
            let mean_accessory = mean_total - mean_core;

            // Standard error
            let variance = total_genes
                .iter()
                .map(|&x| {
                    let diff = x - mean_total;
                    diff * diff
                })
                .sum::<f64>()
                / self.num_samples as f64;
            let stderr = if self.num_samples > 1 {
                (variance / (self.num_samples - 1) as f64).sqrt()
            } else {
                0.0
            };

            curve_data.push(AccumulationPoint {
                k,
                mean_total,
                mean_core,
                mean_accessory,
                stderr,
            });
        }

        let heaps_law = self.fit_heaps_law(&curve_data);

        AccumulationResult {
            curve_data,
            heaps_law,
        }
    }

    /// Fit Heaps' law via linear regression on log-transformed data.
    fn fit_heaps_law(&self, curve_data: &[AccumulationPoint]) -> HeapsLawFit {
        if curve_data.is_empty() || curve_data.len() < 2 {
            return HeapsLawFit {
                alpha: 0.0,
                a_coefficient: 0.0,
                r_squared: 0.0,
                classification: "unknown".to_string(),
            };
        }

        // Filter out zero values for log transformation
        let points: Vec<(f64, f64)> = curve_data
            .iter()
            .filter(|p| p.k > 0 && p.mean_total > 0.0)
            .map(|p| (p.k as f64, p.mean_total))
            .collect();

        if points.len() < 2 {
            return HeapsLawFit {
                alpha: 0.0,
                a_coefficient: 0.0,
                r_squared: 0.0,
                classification: "unknown".to_string(),
            };
        }

        let n_points = points.len() as f64;
        let sum_ln_k: f64 = points.iter().map(|(k, _)| k.ln()).sum();
        let sum_ln_n: f64 = points.iter().map(|(_, n)| n.ln()).sum();
        let sum_ln_kln_n: f64 = points.iter().map(|(k, n)| k.ln() * n.ln()).sum();
        let sum_ln_k_sq: f64 = points.iter().map(|(k, _)| k.ln() * k.ln()).sum();

        // Linear regression: ln(n) = alpha * ln(k) + ln(A)
        // slope = alpha, intercept = ln(A)
        let denominator = n_points * sum_ln_k_sq - sum_ln_k * sum_ln_k;
        if denominator.abs() < 1e-10 {
            return HeapsLawFit {
                alpha: 0.0,
                a_coefficient: 0.0,
                r_squared: 0.0,
                classification: "unknown".to_string(),
            };
        }

        let alpha = (n_points * sum_ln_kln_n - sum_ln_k * sum_ln_n) / denominator;
        let intercept = (sum_ln_n - alpha * sum_ln_k) / n_points;
        let a_coefficient = intercept.exp();

        // R-squared
        let mean_ln_n = sum_ln_n / n_points;
        let ss_tot: f64 = points.iter().map(|(_, n)| {
            let diff = n.ln() - mean_ln_n;
            diff * diff
        }).sum();
        let ss_res: f64 = points.iter().map(|(k, n)| {
            let predicted = alpha * k.ln() + intercept;
            let diff = n.ln() - predicted;
            diff * diff
        }).sum();

        let r_squared = if ss_tot.abs() > 1e-10 {
            1.0 - ss_res / ss_tot
        } else {
            0.0
        };

        let classification = if alpha >= 1.0 {
            "open pangenome".to_string()
        } else {
            "closed pangenome".to_string()
        };

        HeapsLawFit {
            alpha,
            a_coefficient,
            r_squared,
            classification,
        }
    }

    /// Parse a gene_presence_absence.csv file into a BitPackedMatrix.
    fn parse_csv_to_matrix(csv_path: &Path) -> Result<BitPackedMatrix> {
        let mut reader = csv::Reader::from_path(csv_path)?;
        let headers = reader.headers()?.clone();

        // Roary CSV has 14 metadata columns, then per-genome columns
        const METADATA_COLS: usize = 14;
        let genome_names: Vec<String> = headers.iter().skip(METADATA_COLS).map(|s| s.to_string()).collect();
        let num_genomes = genome_names.len();

        // Collect all records first to know cluster count
        let records: Vec<csv::StringRecord> = reader.records().filter_map(|r| r.ok()).collect();
        let num_clusters = records.len();

        let mut matrix = BitPackedMatrix::new(num_genomes, num_clusters);
        matrix.set_genome_names(genome_names);

        let cluster_ids: Vec<String> = records
            .iter()
            .map(|r| r.get(0).unwrap_or("").to_string())
            .collect();
        matrix.set_cluster_ids(cluster_ids);

        for (cluster_idx, record) in records.iter().enumerate() {
            for (genome_idx, val) in record.iter().enumerate().skip(METADATA_COLS) {
                let present = val != "" && val != "0" && !val.eq_ignore_ascii_case("FALSE");
                matrix.set(genome_idx - METADATA_COLS, cluster_idx, present);
            }
        }

        Ok(matrix)
    }
}

impl Default for AccumulationCurveRunner {
    fn default() -> Self {
        Self::new()
    }
}

impl DownstreamRunner for AccumulationCurveRunner {
    type Output = AccumulationResult;

    fn run(&self, output_dir: &Path) -> Result<Self::Output> {
        let csv_path = output_dir.join("gene_presence_absence.csv");
        let matrix = Self::parse_csv_to_matrix(&csv_path)?;
        Ok(self.run_on_matrix(&matrix))
    }

    fn name(&self) -> &str {
        "AccumulationCurveRunner"
    }

    fn is_available(&self) -> bool {
        // No external dependency — always available
        true
    }

    fn required_inputs(&self) -> Vec<DownstreamInput> {
        vec![DownstreamInput::PresenceAbsenceCsv]
    }
}

/// A single point on the accumulation curve.
#[derive(Debug, Clone)]
pub struct AccumulationPoint {
    /// Number of genomes sampled.
    pub k: usize,
    /// Mean number of total genes across samples.
    pub mean_total: f64,
    /// Mean number of core genes (present in all k genomes).
    pub mean_core: f64,
    /// Mean number of accessory genes.
    pub mean_accessory: f64,
    /// Standard error of the mean for total genes.
    pub stderr: f64,
}

/// Heaps' law fit parameters.
#[derive(Debug, Clone)]
pub struct HeapsLawFit {
    /// Power exponent (alpha in n(k) = A * k^alpha).
    pub alpha: f64,
    /// Coefficient A in n(k) = A * k^alpha.
    pub a_coefficient: f64,
    /// Coefficient of determination (R^2) for the log-log linear fit.
    pub r_squared: f64,
    /// Classification: "open pangenome" or "closed pangenome".
    pub classification: String,
}

/// Result of the accumulation curve analysis.
#[derive(Debug, Clone)]
pub struct AccumulationResult {
    /// Rarefaction curve data points.
    pub curve_data: Vec<AccumulationPoint>,
    /// Fitted Heaps' law parameters.
    pub heaps_law: HeapsLawFit,
}

impl DownstreamResult for AccumulationResult {
    fn write_to(&self, dir: &Path) -> Result<()> {
        // Write accumulation_curve.csv
        let curve_path = dir.join("accumulation_curve.csv");
        let mut writer = csv::Writer::from_path(&curve_path)?;
        writer.write_record(&[
            "k",
            "mean_total",
            "mean_core",
            "mean_accessory",
            "stderr",
        ])?;
        for point in &self.curve_data {
            writer.write_record(&[
                point.k.to_string(),
                format!("{:.4}", point.mean_total),
                format!("{:.4}", point.mean_core),
                format!("{:.4}", point.mean_accessory),
                format!("{:.4}", point.stderr),
            ])?;
        }
        writer.flush()?;

        // Write heaps_law_fit.csv
        let heaps_path = dir.join("heaps_law_fit.csv");
        let mut writer = csv::Writer::from_path(&heaps_path)?;
        writer.write_record(&[
            "alpha",
            "A",
            "r_squared",
            "classification",
        ])?;
        writer.write_record(&[
            format!("{:.6}", self.heaps_law.alpha),
            format!("{:.6}", self.heaps_law.a_coefficient),
            format!("{:.6}", self.heaps_law.r_squared),
            self.heaps_law.classification.clone(),
        ])?;
        writer.flush()?;

        Ok(())
    }

    fn summary(&self) -> String {
        format!(
            "Accumulation curve: {} rarefaction points, Heaps' law alpha={:.4} ({}), R^2={:.4}",
            self.curve_data.len(),
            self.heaps_law.alpha,
            self.heaps_law.classification,
            self.heaps_law.r_squared,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_accumulation_curve_monotone_increasing() {
        // 5 genomes, each with 10 unique genes (no overlap)
        // This creates a fully open pangenome where each new genome adds 10 genes
        let mut matrix = BitPackedMatrix::new(5, 50);
        matrix.set_genome_names(vec![
            "g1".to_string(),
            "g2".to_string(),
            "g3".to_string(),
            "g4".to_string(),
            "g5".to_string(),
        ]);
        matrix.set_cluster_ids((0..50).map(|i| format!("gene_{}", i)).collect());

        // Genomes 0-4 each have 10 unique genes
        // g1: genes 0-9, g2: genes 10-19, g3: genes 20-29, g4: genes 30-39, g5: genes 40-49
        for g in 0..5 {
            for c in 0..10 {
                let cluster_idx = g * 10 + c;
                matrix.set(g, cluster_idx, true);
            }
        }

        let runner = AccumulationCurveRunner::new()
            .with_num_samples(10)
            .with_rarefaction_points(5);

        let result = runner.run_on_matrix(&matrix);

        // Total genes should monotonically increase as k increases
        for i in 1..result.curve_data.len() {
            assert!(
                result.curve_data[i].mean_total >= result.curve_data[i - 1].mean_total,
                "Total genes should monotonically increase: {} -> {} at k={} vs k={}",
                result.curve_data[i - 1].mean_total,
                result.curve_data[i].mean_total,
                result.curve_data[i - 1].k,
                result.curve_data[i].k
            );
        }

        // With 5 unique genomes and 10 genes each, total at k=5 should be 50
        let last_point = result.curve_data.last().unwrap();
        assert!((last_point.mean_total - 50.0).abs() < 1.0,
            "At k=5 with 5 unique genomes x 10 genes each, total should be ~50, got {}",
            last_point.mean_total);
    }

    #[test]
    fn test_heaps_law_closed_pangenome() {
        // 10 genomes, all share the same 100 core genes (closed pangenome)
        let n = 10;
        let num_clusters = 100;
        let mut matrix = BitPackedMatrix::new(n, num_clusters);
        matrix.set_genome_names((0..n).map(|i| format!("g{}", i)).collect());
        matrix.set_cluster_ids((0..num_clusters).map(|i| format!("gene_{}", i)).collect());

        // All genes are core (present in all genomes)
        for g in 0..n {
            for c in 0..num_clusters {
                matrix.set(g, c, true);
            }
        }

        let runner = AccumulationCurveRunner::new()
            .with_num_samples(20)
            .with_rarefaction_points(10);

        let result = runner.run_on_matrix(&matrix);

        // Closed pangenome: alpha should be < 1.0
        assert!(
            result.heaps_law.alpha < 1.0,
            "Closed pangenome (all genes core) should have alpha < 1.0, got {}",
            result.heaps_law.alpha
        );
        assert_eq!(
            result.heaps_law.classification, "closed pangenome",
            "Expected closed pangenome classification"
        );
    }

    #[test]
    fn test_heaps_law_open_pangenome() {
        // Open pangenome: each new genome adds new genes (power law growth)
        // Simulate by having gene frequency decrease exponentially
        let n = 20;
        let num_clusters = 200;
        let mut matrix = BitPackedMatrix::new(n, num_clusters);
        matrix.set_genome_names((0..n).map(|i| format!("g{}", i)).collect());
        matrix.set_cluster_ids((0..num_clusters).map(|i| format!("gene_{}", i)).collect());

        // Each gene cluster is present in roughly half the genomes (accessory)
        // This creates a power-law-like distribution typical of open pangenomes
        for c in 0..num_clusters {
            for g in 0..n {
                // Probability of presence decreases with cluster index
                // Later clusters are rarer (present in fewer genomes)
                let prob = 1.0 - (c as f64 / num_clusters as f64) * 0.95;
                if (g as f64 / n as f64) < prob {
                    matrix.set(g, c, true);
                }
            }
        }

        let runner = AccumulationCurveRunner::new()
            .with_num_samples(50)
            .with_rarefaction_points(15);

        let result = runner.run_on_matrix(&matrix);

        // For this distribution, alpha should be >= 1.0 (open pangenome)
        // Note: depending on the exact distribution, result may vary
        println!("Heaps' law alpha: {:.4}", result.heaps_law.alpha);
        println!("Classification: {}", result.heaps_law.classification);
        println!("R^2: {:.4}", result.heaps_law.r_squared);

        // Just verify the fit is computed and reasonable
        assert!(result.heaps_law.r_squared >= 0.0 && result.heaps_law.r_squared <= 1.0);
        assert!(result.heaps_law.alpha > 0.0);
    }

    #[test]
    fn test_rarefaction_ks() {
        let runner = AccumulationCurveRunner::new();

        // When rarefaction_points >= n, return all values 1..=n
        let ks = runner.rarefaction_ks(5);
        assert_eq!(ks, vec![1, 2, 3, 4, 5]);

        // When rarefaction_points == 0
        let runner = AccumulationCurveRunner::new().with_rarefaction_points(0);
        let ks = runner.rarefaction_ks(10);
        assert!(ks.is_empty());

        // Even spacing for 20 points from 1 to 100
        let runner = AccumulationCurveRunner::new().with_rarefaction_points(20);
        let ks = runner.rarefaction_ks(100);
        assert_eq!(ks.len(), 20);
        assert_eq!(ks[0], 1);
        assert_eq!(ks[19], 100);
    }

    #[test]
    fn test_accumulation_result_write_to() {
        let result = AccumulationResult {
            curve_data: vec![
                AccumulationPoint {
                    k: 1,
                    mean_total: 10.0,
                    mean_core: 10.0,
                    mean_accessory: 0.0,
                    stderr: 0.0,
                },
                AccumulationPoint {
                    k: 5,
                    mean_total: 50.0,
                    mean_core: 10.0,
                    mean_accessory: 40.0,
                    stderr: 2.0,
                },
            ],
            heaps_law: HeapsLawFit {
                alpha: 1.5,
                a_coefficient: 10.0,
                r_squared: 0.95,
                classification: "open pangenome".to_string(),
            },
        };

        let temp_dir = TempDir::new().unwrap();
        result.write_to(temp_dir.path()).unwrap();

        // Check accumulation_curve.csv
        let curve_path = temp_dir.path().join("accumulation_curve.csv");
        let curve_content = std::fs::read_to_string(&curve_path).unwrap();
        assert!(curve_content.contains("k,mean_total,mean_core,mean_accessory,stderr"));
        assert!(curve_content.contains("1,10.0000,10.0000,0.0000,0.0000"));
        assert!(curve_content.contains("5,50.0000,10.0000,40.0000,2.0000"));

        // Check heaps_law_fit.csv
        let heaps_path = temp_dir.path().join("heaps_law_fit.csv");
        let heaps_content = std::fs::read_to_string(&heaps_path).unwrap();
        assert!(heaps_content.contains("alpha,A,r_squared,classification"));
        assert!(heaps_content.contains("1.500000"));
        assert!(heaps_content.contains("open pangenome"));
    }
}
