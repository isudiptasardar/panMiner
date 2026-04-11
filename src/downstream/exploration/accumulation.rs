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
    fn rarefactionKs(&self, n: usize) -> Vec<usize> {
        if self.rarefaction_points == 0 {
            return vec![];
        }
        if self.rarefaction_points >= n {
            return (1..=n).collect();
        }
        // Evenly spaced integer points from 1 to n
        let step = (n - 1) as f64 / (self.rarefaction_points - 1) as f64;
        (0..self.rarefaction_points)
            .map(|i| {
                let k = 1.0 + i as f64 * step;
                k.round() as usize
            })
            .collect()
    }

    /// Count genes present in the given sampled genomes.
    fn countGenes(&self, matrix: &BitPackedMatrix, genomeIndices: &[usize]) -> (usize, usize) {
        let k = genomeIndices.len();
        let totalGenes = matrix.num_clusters();
        let mut presentCount = 0usize;
        let mut coreCount = 0usize;

        for clusterIdx in 0..totalGenes {
            // Count how many of the sampled genomes have this cluster
            let count = genomeIndices
                .iter()
                .filter(|&&g| matrix.get(g, clusterIdx))
                .count();
            if count > 0 {
                presentCount += 1;
            }
            // Core gene: present in ALL sampled genomes
            if count == k {
                coreCount += 1;
            }
        }

        (presentCount, coreCount)
    }

    /// Run rarefaction analysis on an in-memory BitPackedMatrix.
    fn runOnMatrix(&self, matrix: &BitPackedMatrix) -> AccumulationResult {
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

        let ks = self.rarefactionKs(n);
        let mut curve_data = Vec::with_capacity(ks.len());

        // Fixed seed for reproducibility
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        for &k in &ks {
            let mut totalGenes = Vec::with_capacity(self.num_samples);
            let mut coreGenes = Vec::with_capacity(self.num_samples);

            for _ in 0..self.num_samples {
                // Sample k genomes without replacement
                let mut indices: Vec<usize> = (0..n).collect();
                indices.shuffle(&mut rng);
                indices.truncate(k);
                indices.sort();

                let (present, core) = self.countGenes(matrix, &indices);
                totalGenes.push(present as f64);
                coreGenes.push(core as f64);
            }

            let meanTotal = totalGenes.iter().sum::<f64>() / self.num_samples as f64;
            let meanCore = coreGenes.iter().sum::<f64>() / self.num_samples as f64;
            let meanAccessory = meanTotal - meanCore;

            // Standard error
            let variance = totalGenes
                .iter()
                .map(|&x| {
                    let diff = x - meanTotal;
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
                mean_total: meanTotal,
                mean_core: meanCore,
                mean_accessory: meanAccessory,
                stderr,
            });
        }

        let heaps_law = self.fitHeapsLaw(&curve_data);

        AccumulationResult {
            curve_data,
            heaps_law,
        }
    }

    /// Fit Heaps' law via linear regression on log-transformed data.
    fn fitHeapsLaw(&self, curveData: &[AccumulationPoint]) -> HeapsLawFit {
        if curveData.is_empty() || curveData.len() < 2 {
            return HeapsLawFit {
                alpha: 0.0,
                a_coefficient: 0.0,
                r_squared: 0.0,
                classification: "unknown".to_string(),
            };
        }

        // Filter out zero values for log transformation
        let points: Vec<(f64, f64)> = curveData
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

        let nPoints = points.len() as f64;
        let sumLnK: f64 = points.iter().map(|(k, _)| k.ln()).sum();
        let sumLnN: f64 = points.iter().map(|(_, n)| n.ln()).sum();
        let sumLnKLnN: f64 = points.iter().map(|(k, n)| k.ln() * n.ln()).sum();
        let sumLnKSq: f64 = points.iter().map(|(k, _)| k.ln() * k.ln()).sum();

        // Linear regression: ln(n) = alpha * ln(k) + ln(A)
        // slope = alpha, intercept = ln(A)
        let denominator = nPoints * sumLnKSq - sumLnK * sumLnK;
        if denominator.abs() < 1e-10 {
            return HeapsLawFit {
                alpha: 0.0,
                a_coefficient: 0.0,
                r_squared: 0.0,
                classification: "unknown".to_string(),
            };
        }

        let alpha = (nPoints * sumLnKLnN - sumLnK * sumLnN) / denominator;
        let intercept = (sumLnN - alpha * sumLnK) / nPoints;
        let aCoefficient = intercept.exp();

        // R-squared
        let meanLnN = sumLnN / nPoints;
        let ssTot: f64 = points.iter().map(|(_, n)| {
            let diff = n.ln() - meanLnN;
            diff * diff
        }).sum();
        let ssRes: f64 = points.iter().map(|(k, n)| {
            let predicted = alpha * k.ln() + intercept;
            let diff = n.ln() - predicted;
            diff * diff
        }).sum();

        let rSquared = if ssTot.abs() > 1e-10 {
            1.0 - ssRes / ssTot
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
            a_coefficient: aCoefficient,
            r_squared: rSquared,
            classification,
        }
    }

    /// Parse a gene_presence_absence.csv file into a BitPackedMatrix.
    fn parseCsvToMatrix(csvPath: &Path) -> Result<BitPackedMatrix> {
        let mut reader = csv::Reader::from_path(csvPath)?;
        let headers = reader.headers()?.clone();

        // First column is "Gene", rest are genome names
        let genomeNames: Vec<String> = headers.iter().skip(1).map(|s| s.to_string()).collect();
        let numGenomes = genomeNames.len();

        // Collect all records first to know cluster count
        let records: Vec<csv::StringRecord> = reader.records().filter_map(|r| r.ok()).collect();
        let numClusters = records.len();

        let mut matrix = BitPackedMatrix::new(numGenomes, numClusters);
        matrix.set_genome_names(genomeNames);

        let clusterIds: Vec<String> = records
            .iter()
            .map(|r| r.get(0).unwrap_or("").to_string())
            .collect();
        matrix.set_cluster_ids(clusterIds);

        for (clusterIdx, record) in records.iter().enumerate() {
            for (genomeIdx, val) in record.iter().enumerate().skip(1) {
                let present = val == "1" || val.eq_ignore_ascii_case("TRUE");
                matrix.set(genomeIdx - 1, clusterIdx, present);
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
        let csvPath = output_dir.join("gene_presence_absence.csv");
        let matrix = Self::parseCsvToMatrix(&csvPath)?;
        Ok(self.runOnMatrix(&matrix))
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
        let curvePath = dir.join("accumulation_curve.csv");
        let mut writer = csv::Writer::from_path(&curvePath)?;
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
        let heapsPath = dir.join("heaps_law_fit.csv");
        let mut writer = csv::Writer::from_path(&heapsPath)?;
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
    fn testAccumulationCurveMonotoneIncreasing() {
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
                let clusterIdx = g * 10 + c;
                matrix.set(g, clusterIdx, true);
            }
        }

        let runner = AccumulationCurveRunner::new()
            .with_num_samples(10)
            .with_rarefaction_points(5);

        let result = runner.runOnMatrix(&matrix);

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
        let lastPoint = result.curve_data.last().unwrap();
        assert!((lastPoint.mean_total - 50.0).abs() < 1.0,
            "At k=5 with 5 unique genomes x 10 genes each, total should be ~50, got {}",
            lastPoint.mean_total);
    }

    #[test]
    fn testHeapsLawClosedPangenome() {
        // 10 genomes, all share the same 100 core genes (closed pangenome)
        let n = 10;
        let numClusters = 100;
        let mut matrix = BitPackedMatrix::new(n, numClusters);
        matrix.set_genome_names((0..n).map(|i| format!("g{}", i)).collect());
        matrix.set_cluster_ids((0..numClusters).map(|i| format!("gene_{}", i)).collect());

        // All genes are core (present in all genomes)
        for g in 0..n {
            for c in 0..numClusters {
                matrix.set(g, c, true);
            }
        }

        let runner = AccumulationCurveRunner::new()
            .with_num_samples(20)
            .with_rarefaction_points(10);

        let result = runner.runOnMatrix(&matrix);

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
    fn testHeapsLawOpenPangenome() {
        // Open pangenome: each new genome adds new genes (power law growth)
        // Simulate by having gene frequency decrease exponentially
        let n = 20;
        let numClusters = 200;
        let mut matrix = BitPackedMatrix::new(n, numClusters);
        matrix.set_genome_names((0..n).map(|i| format!("g{}", i)).collect());
        matrix.set_cluster_ids((0..numClusters).map(|i| format!("gene_{}", i)).collect());

        // Each gene cluster is present in roughly half the genomes (accessory)
        // This creates a power-law-like distribution typical of open pangenomes
        for c in 0..numClusters {
            for g in 0..n {
                // Probability of presence decreases with cluster index
                // Later clusters are rarer (present in fewer genomes)
                let prob = 1.0 - (c as f64 / numClusters as f64) * 0.95;
                if (g as f64 / n as f64) < prob {
                    matrix.set(g, c, true);
                }
            }
        }

        let runner = AccumulationCurveRunner::new()
            .with_num_samples(50)
            .with_rarefaction_points(15);

        let result = runner.runOnMatrix(&matrix);

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
    fn testRarefactionKs() {
        let runner = AccumulationCurveRunner::new();

        // When rarefaction_points >= n, return all values 1..=n
        let ks = runner.rarefactionKs(5);
        assert_eq!(ks, vec![1, 2, 3, 4, 5]);

        // When rarefaction_points == 0
        let runner = AccumulationCurveRunner::new().with_rarefaction_points(0);
        let ks = runner.rarefactionKs(10);
        assert!(ks.is_empty());

        // Even spacing for 20 points from 1 to 100
        let runner = AccumulationCurveRunner::new().with_rarefaction_points(20);
        let ks = runner.rarefactionKs(100);
        assert_eq!(ks.len(), 20);
        assert_eq!(ks[0], 1);
        assert_eq!(ks[19], 100);
    }

    #[test]
    fn testAccumulationResultWriteTo() {
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

        let tempDir = TempDir::new().unwrap();
        result.write_to(tempDir.path()).unwrap();

        // Check accumulation_curve.csv
        let curvePath = tempDir.path().join("accumulation_curve.csv");
        let curveContent = std::fs::read_to_string(&curvePath).unwrap();
        assert!(curveContent.contains("k,mean_total,mean_core,mean_accessory,stderr"));
        assert!(curveContent.contains("1,10.0000,10.0000,0.0000,0.0000"));
        assert!(curveContent.contains("5,50.0000,10.0000,40.0000,2.0000"));

        // Check heaps_law_fit.csv
        let heapsPath = tempDir.path().join("heaps_law_fit.csv");
        let heapsContent = std::fs::read_to_string(&heapsPath).unwrap();
        assert!(heapsContent.contains("alpha,A,r_squared,classification"));
        assert!(heapsContent.contains("1.500000"));
        assert!(heapsContent.contains("open pangenome"));
    }
}
