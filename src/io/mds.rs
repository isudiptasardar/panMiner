//! Classical MDS (metric multidimensional scaling) projection.
//!
//! This module provides pure-Rust MDS projection from a distance matrix into 2D
//! space, used for QC visualization scatter plots. No external dependencies are
//! required — the implementation uses power iteration for eigendecomposition.

use std::io::Write;

use crate::error::Result;

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

/// Compute MDS projection from a distance matrix with provided labels.
///
/// Uses classical MDS (metric MDS) to project genomes into 2D space
/// for visualization. Power iteration is used for eigendecomposition.
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
}