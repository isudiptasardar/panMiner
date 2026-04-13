//! Sourmash-based genome distance estimation and MDS projection.
//!
//! This module provides:
//! - `SourmashRunner`: Subprocess runner for computing MinHash distance matrices
//!   via the `sourmash` CLI (always compiled, requires sourmash installed)
//! - `MdsProjection`: Classical MDS projection from a distance matrix into 2D
//! - `compute_mds` / `compute_mds_with_labels`: MDS computation (pure Rust, always available)
//! - `compute_distance_matrix`: Feature-gated direct sourmash API integration
//!
//! The core MDS projection and `SourmashRunner` are always available.
//! The `sourmash` feature gate only controls the direct library integration
//! (`compute_distance_matrix`), which requires the sourmash Rust bindings.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::Result;

/// Sourmash subprocess runner for computing MinHash distance matrices.
///
/// Uses `sourmash sketch` to create signatures and `sourmash compare`
/// to compute pairwise distances between genomes.
pub struct SourmashRunner {
    /// Path to sourmash binary
    sourmash_path: PathBuf,
}

impl SourmashRunner {
    /// Create a new SourmashRunner with an explicit path.
    pub fn new(sourmash_path: PathBuf) -> Self {
        Self { sourmash_path }
    }

    /// Detect if sourmash is installed on the system.
    ///
    /// Returns `Some(SourmashRunner)` if `sourmash --version` succeeds,
    /// `None` otherwise.
    pub fn detect() -> Option<Self> {
        let path = which_sourmash()?;
        Some(Self { sourmash_path: path })
    }

    /// Get the path to the sourmash binary.
    pub fn path(&self) -> &Path {
        &self.sourmash_path
    }

    /// Compute a pairwise distance matrix for all genomes using sourmash.
    ///
    /// Workflow:
    /// 1. Create MinHash sketches for each genome using `sourmash sketch`
    /// 2. Compute pairwise distances using `sourmash compare`
    /// 3. Parse the output into a symmetric distance matrix
    ///
    /// Returns a matrix where result[i][j] is the distance between genome i
    /// and genome j (0.0 = identical, 1.0 = completely different).
    pub fn compute_distance_matrix(&self, genomes: &[PathBuf]) -> Result<Vec<Vec<f64>>> {
        if genomes.is_empty() {
            return Ok(Vec::new());
        }

        let temp_dir = tempfile::tempdir()?;
        let sig_dir = temp_dir.path().join("signatures");
        std::fs::create_dir_all(&sig_dir)?;

        // Step 1: Create sketches for each genome
        tracing::info!("Creating MinHash sketches for {} genomes", genomes.len());
        for genome_path in genomes {
            let genome_name = genome_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown");
            let sig_path = sig_dir.join(format!("{}.sig.gz", genome_name));

            let output = Command::new(&self.sourmash_path)
                .arg("sketch")
                .arg("dna")
                .arg("-p")
                .arg("1") // single thread for sketching
                .arg("-o")
                .arg(&sig_path)
                .arg(genome_path)
                .output()
                .map_err(|e| crate::Error::ExternalTool(format!("sourmash sketch failed: {}", e)))?;

            if !output.status.success() {
                let stderr = String::from_utf8_lossy(&output.stderr);
                return Err(crate::Error::ExternalTool(format!(
                    "sourmash sketch failed for {}: {}",
                    genome_name,
                    stderr.trim()
                )));
            }
        }

        // Step 2: Collect all signature files
        let sig_files: Vec<PathBuf> = std::fs::read_dir(&sig_dir)?
            .filter_map(|entry| {
                let e = entry.ok()?;
                let path = e.path();
                if path.extension().map_or(false, |ext| ext == "gz" || ext == "sig") {
                    Some(path)
                } else {
                    None
                }
            })
            .collect();

        if sig_files.len() != genomes.len() {
            return Err(crate::Error::Output(format!(
                "Expected {} signature files, found {}",
                genomes.len(),
                sig_files.len()
            )));
        }

        // Step 3: Run sourmash compare
        let compare_output_path = temp_dir.path().join("compare.csv");
        let sig_list_path = temp_dir.path().join("sig_list.txt");
        {
            let mut list_file = std::fs::File::create(&sig_list_path)?;
            for sig in &sig_files {
                writeln!(list_file, "{}", sig.display())?;
            }
        }

        tracing::info!("Computing pairwise distances for {} genomes", genomes.len());
        let output = Command::new(&self.sourmash_path)
            .arg("compare")
            .arg("-k")
            .arg("31")
            .arg("--dna")
            .arg("-o")
            .arg(&compare_output_path)
            .arg("--from-file")
            .arg(&sig_list_path)
            .output()
            .map_err(|e| crate::Error::ExternalTool(format!("sourmash compare failed: {}", e)))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(crate::Error::ExternalTool(format!(
                "sourmash compare failed: {}",
                stderr.trim()
            )));
        }

        // Step 4: Parse the compare output (CSV with header of genome names)
        parse_sourmash_compare(&compare_output_path, genomes)
    }
}

/// Find the sourmash binary on PATH.
fn which_sourmash() -> Option<PathBuf> {
    which::which("sourmash").ok()
}

/// Parse sourmash compare output into a distance matrix.
///
/// The sourmash compare output CSV has format:
/// ```text
/// ,genome_0,genome_1,...
/// genome_0,0.00,0.05,...
/// genome_1,0.05,0.00,...
/// ...
/// ```
///
/// The values are Jaccard similarities (0=identical, 1=completely different).
fn parse_sourmash_compare(path: &Path, genomes: &[PathBuf]) -> Result<Vec<Vec<f64>>> {
    let content = std::fs::read_to_string(path)?;
    let mut lines = content.lines();

    // Skip header line
    let _header = lines.next().ok_or_else(|| {
        crate::Error::Output("sourmash compare output is empty".to_string())
    })?;

    let n = genomes.len();
    let mut matrix = vec![vec![0.0; n]; n];

    for (i, line) in lines.enumerate() {
        if i >= n {
            break;
        }
        let parts: Vec<&str> = line.split(',').collect();
        // First column is the row label, then distance values
        for (j, val_str) in parts.iter().skip(1).enumerate() {
            if j < n {
                if let Ok(similarity) = val_str.trim().parse::<f64>() {
                    // sourmash compare outputs similarity (1-distance),
                    // we convert to distance
                    let distance = 1.0 - similarity;
                    matrix[i][j] = distance;
                    matrix[j][i] = distance;
                }
            }
        }
        // Diagonal is 0 (self-distance)
        matrix[i][i] = 0.0;
    }

    Ok(matrix)
}

/// Find the length from the header (number of columns minus label column).
#[allow(dead_code)]
fn parse_header_count(header: &str) -> usize {
    header.split(',').count().saturating_sub(1)
}

/// Result of an MDS projection.
#[derive(Debug, Clone)]
pub struct MdsProjection {
    /// 2D coordinates for each genome.
    pub coordinates: Vec<(f64, f64)>,
    /// Genome labels.
    pub labels: Vec<String>,
}

impl MdsProjection {
    /// Write MDS coordinates to a CSV file.
    pub fn write_csv(&self, path: &std::path::Path) -> Result<()> {
        let file = std::fs::File::create(path)?;
        let mut writer = std::io::BufWriter::new(file);
        writeln!(writer, "genome,x,y")?;
        for (i, label) in self.labels.iter().enumerate() {
            let (x, y) = self.coordinates[i];
            writeln!(writer, "{},{:.6},{:.6}", label, x, y)?;
        }
        Ok(())
    }

    /// Generate a scatter plot HTML using d3.js.
    pub fn write_html(&self, path: &std::path::Path) -> Result<()> {
        let file = std::fs::File::create(path)?;
        let mut writer = std::io::BufWriter::new(file);
        let data: String = serde_json::to_string(&self.coordinates).unwrap_or_else(|_| "[]".to_string());
        let labels_json: String = serde_json::to_string(&self.labels).unwrap_or_else(|_| "[]".to_string());

        writeln!(writer, "<!DOCTYPE html>")?;
        writeln!(writer, "<html lang='en'>")?;
        writeln!(writer, "<head>")?;
        writeln!(writer, "  <meta charset='UTF-8'>")?;
        writeln!(writer, "  <meta name='viewport' content='width=device-width, initial-scale=1.0'>")?;
        writeln!(writer, "  <title>PanMiner MDS Projection</title>")?;
        writeln!(writer, "  <script src='https://d3js.org/d3.v7.min.js'></script>")?;
        writeln!(writer, "  <style>")?;
        writeln!(writer, "    body {{ font-family: sans-serif; background: #1a1a2e; color: #eee; margin: 0; padding: 20px; }}")?;
        writeln!(writer, "    h1 {{ color: #e94560; margin-bottom: 20px; }}")?;
        writeln!(writer, "    .chart {{ background: #0a0a15; border-radius: 8px; padding: 20px; }}")?;
        writeln!(writer, "    .point {{ fill: #e94560; cursor: pointer; }}")?;
        writeln!(writer, "    .point:hover {{ fill: #4ecdc4; }}")?;
        writeln!(writer, "    .axis-label {{ fill: #888; font-size: 12px; }}")?;
        writeln!(writer, "    .tooltip {{ position: absolute; background: rgba(0,0,0,0.9); padding: 8px 12px; border-radius: 4px; font-size: 13px; pointer-events: none; }}")?;
        writeln!(writer, "  </style>")?;
        writeln!(writer, "</head>")?;
        writeln!(writer, "<body>")?;
        writeln!(writer, "  <h1>PanMiner MDS Projection</h1>")?;
        writeln!(writer, "  <div class='chart' id='chart'></div>")?;
        writeln!(writer, "  <div class='tooltip' id='tooltip' style='display:none;'></div>")?;
        writeln!(writer, "  <script>")?;
        writeln!(writer, "    const coords = {};", data)?;
        writeln!(writer, "    const labels = {};", labels_json)?;
        writeln!(writer, "    const tooltip = document.getElementById('tooltip');")?;
        writeln!(writer, "    const chart = document.getElementById('chart');")?;
        writeln!(writer, "    const w = chart.clientWidth || 800;")?;
        writeln!(writer, "    const h = 500;")?;
        writeln!(writer, "    const margin = {{ top: 20, right: 20, bottom: 40, left: 50 }};")?;
        writeln!(writer, "    const svg = d3.select(chart).append('svg').attr('width', w).attr('height', h);")?;
        writeln!(writer, "    const xExt = d3.extent(coords, d => d[0]);")?;
        writeln!(writer, "    const yExt = d3.extent(coords, d => d[1]);")?;
        writeln!(writer, "    const xScale = d3.scaleLinear().domain(xExt).range([margin.left, w - margin.right]);")?;
        writeln!(writer, "    const yScale = d3.scaleLinear().domain(yExt).range([h - margin.bottom, margin.top]);")?;
        writeln!(writer, "    svg.selectAll('.point').data(coords).enter().append('circle')")?;
        writeln!(writer, "      .attr('class', 'point').attr('cx', d => xScale(d[0])).attr('cy', d => yScale(d[1]))")?;
        writeln!(writer, "      .attr('r', 6).on('mouseover', (event, d) => {{")?;
        writeln!(writer, "        const i = coords.indexOf(d);")?;
        writeln!(writer, "        tooltip.style.display = 'block';")?;
        writeln!(writer, "        tooltip.innerHTML = labels[i];")?;
        writeln!(writer, "        tooltip.style.left = (event.pageX + 10) + 'px';")?;
        writeln!(writer, "        tooltip.style.top = (event.pageY - 10) + 'px';")?;
        writeln!(writer, "      }}).on('mouseout', () => {{ tooltip.style.display = 'none'; }});")?;
        writeln!(writer, "    svg.append('text').attr('class', 'axis-label').attr('x', w/2).attr('y', h - 5)")?;
        writeln!(writer, "      .attr('text-anchor', 'middle').text('MDS Dimension 1');")?;
        writeln!(writer, "    svg.append('text').attr('class', 'axis-label').attr('transform', 'rotate(-90)')")?;
        writeln!(writer, "      .attr('x', -h/2).attr('y', 15).attr('text-anchor', 'middle').text('MDS Dimension 2');")?;
        writeln!(writer, "  </script>")?;
        writeln!(writer, "</body>")?;
        writeln!(writer, "</html>")?;
        Ok(())
    }
}

/// Compute a pairwise distance matrix using sourmash.
///
/// When the `sourmash` feature is not enabled, this logs a warning
/// and returns an error.
#[allow(dead_code)]
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

/// Compute MDS projection from a distance matrix with provided labels.
///
/// Uses classical MDS (metric MDS) to project genomes into 2D space
/// for visualization.
pub fn compute_mds_with_labels(distances: &[Vec<f64>], labels: &[String]) -> Result<MdsProjection> {
    let n = distances.len();
    if n == 0 {
        return Ok(MdsProjection {
            coordinates: Vec::new(),
            labels: Vec::new(),
        });
    }

    if labels.len() != n {
        return Err(crate::Error::Output(format!(
            "Labels count {} does not match distance matrix size {}",
            labels.len(),
            n
        )));
    }

    // Classical MDS: eigen-decomposition of centered distance matrix

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
    for _dim in 0..2 {
        let mut v = vec![1.0 / (n as f64).sqrt(); n];
        for _ in 0..100 {
            let mut v_new = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    v_new[i] += b[i][j] * v[j];
                }
            }
            let norm: f64 = v_new.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm > 0.0 {
                for x in v_new.iter_mut() {
                    *x /= norm;
                }
            }
            v = v_new;
        }
        let mut bv = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                bv[i] += b[i][j] * v[j];
            }
        }
        let eigenvalue: f64 = v.iter().zip(bv.iter()).map(|(vi, bvi)| vi * bvi).sum();
        let sqrt_eigenvalue = if eigenvalue > 0.0 { eigenvalue.sqrt() } else { 0.0 };

        for i in 0..n {
            for j in 0..n {
                b[i][j] -= eigenvalue * v[i] * v[j];
            }
        }

        coords.push((v, sqrt_eigenvalue));
    }

    let mut coordinates = Vec::new();
    for i in 0..n {
        let x = coords[0].0[i] * coords[0].1;
        let y = coords[1].0[i] * coords[1].1;
        coordinates.push((x, y));
    }

    Ok(MdsProjection {
        coordinates,
        labels: labels.to_vec(),
    })
}

/// Compute MDS projection from a distance matrix (placeholder labels).
///
/// NOTE: Prefer `compute_mds_with_labels` to provide real genome labels.
pub fn compute_mds(distances: &[Vec<f64>]) -> Result<MdsProjection> {
    let n = distances.len();
    if n == 0 {
        return Ok(MdsProjection {
            coordinates: Vec::new(),
            labels: Vec::new(),
        });
    }

    let labels: Vec<String> = (0..n).map(|i| format!("genome_{}", i)).collect();
    compute_mds_with_labels(distances, &labels)
}

#[cfg(feature = "sourmash")]
fn compute_distance_matrix_sourmash(genomes: &[PathBuf]) -> Result<Vec<Vec<f64>>> {
    let runner = SourmashRunner::detect().ok_or_else(|| {
        crate::Error::ExternalTool(
            "sourmash not found. Install it with: pip install sourmash".to_string()
        )
    })?;
    runner.compute_distance_matrix(genomes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_mds_with_labels_empty() {
        let result = compute_mds_with_labels(&[], &[]).unwrap();
        assert!(result.coordinates.is_empty());
        assert!(result.labels.is_empty());
    }

    #[test]
    fn test_compute_mds_with_labels_identity() {
        let distances = vec![
            vec![0.0, 0.05, 0.1],
            vec![0.05, 0.0, 0.08],
            vec![0.1, 0.08, 0.0],
        ];
        let labels = vec!["genome_a".to_string(), "genome_b".to_string(), "genome_c".to_string()];
        let result = compute_mds_with_labels(&distances, &labels).unwrap();
        assert_eq!(result.coordinates.len(), 3);
        assert_eq!(result.labels, labels);
    }

    #[test]
    fn test_compute_mds_with_labels_mismatch() {
        let distances = vec![vec![0.0, 0.05], vec![0.05, 0.0]];
        let labels = vec!["genome_a".to_string()];
        let result = compute_mds_with_labels(&distances, &labels);
        assert!(result.is_err());
    }

    #[test]
    fn test_mds_projection_write_csv() {
        let mds = MdsProjection {
            coordinates: vec![(0.1, 0.2), (0.3, 0.4)],
            labels: vec!["g1".to_string(), "g2".to_string()],
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("mds.csv");
        mds.write_csv(&path).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("genome,x,y"));
        assert!(content.contains("g1,0.100000,0.200000"));
    }

    #[test]
    fn test_mds_projection_write_html() {
        let mds = MdsProjection {
            coordinates: vec![(0.1, 0.2), (0.3, 0.4)],
            labels: vec!["g1".to_string(), "g2".to_string()],
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("mds.html");
        mds.write_html(&path).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("PanMiner MDS Projection"));
        assert!(content.contains("d3.v7"));
    }

    #[test]
    fn test_sourmash_runner_creation() {
        let runner = SourmashRunner::new(PathBuf::from("/usr/bin/sourmash"));
        assert_eq!(runner.path(), std::path::Path::new("/usr/bin/sourmash"));
    }

    #[test]
    fn test_sourmash_detect() {
        // This just tests that detect() doesn't panic
        let _ = SourmashRunner::detect();
    }

    #[test]
    fn test_parse_sourmash_compare() {
        let temp_dir = tempfile::tempdir().unwrap();
        let csv_path = temp_dir.path().join("compare.csv");
        let csv_content = ",genome_0,genome_1,genome_2\n\
                           genome_0,1.0,0.95,0.88\n\
                           genome_1,0.95,1.0,0.91\n\
                           genome_2,0.88,0.91,1.0\n";
        std::fs::write(&csv_path, csv_content).unwrap();

        let genomes: Vec<PathBuf> = vec![
            PathBuf::from("genome_0.fasta"),
            PathBuf::from("genome_1.fasta"),
            PathBuf::from("genome_2.fasta"),
        ];

        let matrix = parse_sourmash_compare(&csv_path, &genomes).unwrap();

        // Diagonal should be 0
        assert!((matrix[0][0]).abs() < 0.001);
        assert!((matrix[1][1]).abs() < 0.001);
        assert!((matrix[2][2]).abs() < 0.001);

        // genome_0 vs genome_1: distance = 1 - 0.95 = 0.05
        assert!((matrix[0][1] - 0.05).abs() < 0.001);
        assert!((matrix[1][0] - 0.05).abs() < 0.001);

        // genome_0 vs genome_2: distance = 1 - 0.88 = 0.12
        assert!((matrix[0][2] - 0.12).abs() < 0.001);
    }

    #[test]
    fn test_compute_distance_matrix_no_sourmash() {
        // When sourmash feature is not enabled, should return an error
        #[cfg(not(feature = "sourmash"))]
        {
            let result = compute_distance_matrix(&[]);
            assert!(result.is_err());
        }
    }
}
