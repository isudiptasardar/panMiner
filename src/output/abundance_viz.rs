//! Gene abundance visualization module.
//!
//! Generates an interactive HTML report with:
//! - U-shape gene frequency bar chart (x: #genomes, y: #gene families)
//! - Rarefaction curve (x: #genomes added, y: cumulative pangenome size)
//! - Core/soft-core/shell/cloud partition bars
//!
//! Follows the same D3.js HTML pattern as `qc_viz.rs`.

use crate::error::Result;
use std::io::Write;
use std::path::Path;

/// Partition counts for gene families by prevalence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PartitionCounts {
    /// Gene families present in >= 99% of genomes
    pub core: usize,
    /// Gene families present in 95-99% of genomes
    pub soft_core: usize,
    /// Gene families present in 15-95% of genomes
    pub shell: usize,
    /// Gene families present in < 15% of genomes
    pub cloud: usize,
}

/// Writer for gene abundance HTML reports.
pub struct AbundanceVizWriter;

impl AbundanceVizWriter {
    /// Compute partition counts from presence/absence histogram.
    ///
    /// `presence_counts` is a histogram where each element is
    /// `(num_genomes_present, num_gene_families_with_that_count)`.
    /// `total_genomes` is the total number of genomes in the analysis.
    pub fn compute_partitions(
        presence_counts: &[(usize, usize)],
        total_genomes: usize,
    ) -> PartitionCounts {
        if total_genomes == 0 {
            return PartitionCounts {
                core: 0,
                soft_core: 0,
                shell: 0,
                cloud: 0,
            };
        }

        let mut core = 0usize;
        let mut soft_core = 0usize;
        let mut shell = 0usize;
        let mut cloud = 0usize;

        for &(n_genomes, n_families) in presence_counts {
            let fraction = n_genomes as f64 / total_genomes as f64;
            if fraction >= 0.99 {
                core += n_families;
            } else if fraction >= 0.95 {
                soft_core += n_families;
            } else if fraction >= 0.15 {
                shell += n_families;
            } else {
                cloud += n_families;
            }
        }

        PartitionCounts {
            core,
            soft_core,
            shell,
            cloud,
        }
    }

    /// Generate HTML report with gene frequency and rarefaction plots.
    ///
    /// - `presence_counts`: histogram of (num_genomes, num_gene_families)
    /// - `rarefaction_data`: (n_genomes_added, cumulative_pangenome_size)
    /// - `total_genomes`: total number of genomes
    /// - `total_clusters`: total number of gene clusters
    /// - `output_path`: path to write the HTML file
    pub fn write_report(
        presence_counts: &[(usize, usize)],
        rarefaction_data: &[(usize, usize)],
        total_genomes: usize,
        total_clusters: usize,
        output_path: &Path,
    ) -> Result<()> {
        let partitions = Self::compute_partitions(presence_counts, total_genomes);

        let freq_data_json = serde_json::to_string(presence_counts)
            .unwrap_or_else(|_| "[]".to_string());
        let rarefaction_json = serde_json::to_string(rarefaction_data)
            .unwrap_or_else(|_| "[]".to_string());

        let file = std::fs::File::create(output_path)?;
        let mut writer = std::io::BufWriter::new(file);

        writeln!(writer, "<!DOCTYPE html>")?;
        writeln!(writer, "<html lang='en'>")?;
        writeln!(writer, "<head>")?;
        writeln!(writer, "  <meta charset='UTF-8'>")?;
        writeln!(writer, "  <meta name='viewport' content='width=device-width, initial-scale=1.0'>")?;
        writeln!(writer, "  <title>PanMiner Gene Abundance Report</title>")?;
        writeln!(writer, "  <script src='https://d3js.org/d3.v7.min.js'></script>")?;
        writeln!(writer, "  <style>")?;
        writeln!(writer, "    * {{ margin: 0; padding: 0; box-sizing: border-box; }}")?;
        writeln!(writer, "    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; background: #1a1a2e; color: #eee; padding: 20px; }}")?;
        writeln!(writer, "    .header {{ background: rgba(255,255,255,0.05); padding: 20px 30px; border-radius: 8px; margin-bottom: 20px; }}")?;
        writeln!(writer, "    .header h1 {{ color: #e94560; font-size: 24px; margin-bottom: 8px; }}")?;
        writeln!(writer, "    .stats {{ display: flex; gap: 30px; flex-wrap: wrap; margin-top: 10px; }}")?;
        writeln!(writer, "    .stat-item {{ display: flex; align-items: center; gap: 8px; }}")?;
        writeln!(writer, "    .stat-value {{ font-size: 20px; font-weight: 700; color: #4ecdc4; }}")?;
        writeln!(writer, "    .stat-label {{ color: #888; font-size: 14px; }}")?;
        writeln!(writer, "    .section {{ background: rgba(255,255,255,0.03); border-radius: 8px; padding: 20px; margin-bottom: 20px; }}")?;
        writeln!(writer, "    .section h2 {{ color: #e94560; font-size: 18px; margin-bottom: 15px; border-bottom: 1px solid rgba(255,255,255,0.1); padding-bottom: 8px; }}")?;
        writeln!(writer, "    .chart {{ background: #0a0a15; border-radius: 6px; padding: 15px; }}")?;
        writeln!(writer, "    .axis-label {{ fill: #666; font-size: 12px; }}")?;
        writeln!(writer, "    .bar-freq {{ fill: #4ecdc4; }}")?;
        writeln!(writer, "    .bar-freq:hover {{ fill: #e94560; }}")?;
        writeln!(writer, "    .rarefaction-line {{ fill: none; stroke: #4ecdc4; stroke-width: 2.5; }}")?;
        writeln!(writer, "    .rarefaction-point {{ fill: #e94560; }}")?;
        writeln!(writer, "    .partitions {{ display: flex; gap: 12px; flex-wrap: wrap; margin-top: 10px; }}")?;
        writeln!(writer, "    .partition {{ flex: 1; min-width: 120px; border-radius: 6px; padding: 14px 18px; text-align: center; }}")?;
        writeln!(writer, "    .partition h3 {{ font-size: 13px; text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 6px; }}")?;
        writeln!(writer, "    .partition .count {{ font-size: 28px; font-weight: 700; }}")?;
        writeln!(writer, "    .p-core {{ background: rgba(78,205,196,0.15); border: 1px solid rgba(78,205,196,0.3); }}")?;
        writeln!(writer, "    .p-core h3 {{ color: #4ecdc4; }}")?;
        writeln!(writer, "    .p-core .count {{ color: #4ecdc4; }}")?;
        writeln!(writer, "    .p-soft {{ background: rgba(255,165,0,0.15); border: 1px solid rgba(255,165,0,0.3); }}")?;
        writeln!(writer, "    .p-soft h3 {{ color: #ffa500; }}")?;
        writeln!(writer, "    .p-soft .count {{ color: #ffa500; }}")?;
        writeln!(writer, "    .p-shell {{ background: rgba(233,69,96,0.15); border: 1px solid rgba(233,69,96,0.3); }}")?;
        writeln!(writer, "    .p-shell h3 {{ color: #e94560; }}")?;
        writeln!(writer, "    .p-shell .count {{ color: #e94560; }}")?;
        writeln!(writer, "    .p-cloud {{ background: rgba(150,150,150,0.15); border: 1px solid rgba(150,150,150,0.3); }}")?;
        writeln!(writer, "    .p-cloud h3 {{ color: #999; }}")?;
        writeln!(writer, "    .p-cloud .count {{ color: #999; }}")?;
        writeln!(writer, "  </style>")?;
        writeln!(writer, "</head>")?;
        writeln!(writer, "<body>")?;

        // Header with summary stats
        writeln!(writer, "  <div class='header'>")?;
        writeln!(writer, "    <h1>PanMiner Gene Abundance Report</h1>")?;
        writeln!(writer, "    <div class='stats'>")?;
        writeln!(writer, "      <div class='stat-item'><div class='stat-value'>{}</div><div class='stat-label'>Total Genomes</div></div>", total_genomes)?;
        writeln!(writer, "      <div class='stat-item'><div class='stat-value'>{}</div><div class='stat-label'>Gene Families</div></div>", total_clusters)?;
        writeln!(writer, "    </div>")?;
        writeln!(writer, "  </div>")?;

        // Partition bars
        writeln!(writer, "  <div class='section'>")?;
        writeln!(writer, "    <h2>Gene Family Partition</h2>")?;
        writeln!(writer, "    <div class='partitions'>")?;
        writeln!(writer, "      <div class='partition p-core'><h3>Core (&ge;99%)</h3><div class='count'>{}</div></div>", partitions.core)?;
        writeln!(writer, "      <div class='partition p-soft'><h3>Soft-core (95-99%)</h3><div class='count'>{}</div></div>", partitions.soft_core)?;
        writeln!(writer, "      <div class='partition p-shell'><h3>Shell (15-95%)</h3><div class='count'>{}</div></div>", partitions.shell)?;
        writeln!(writer, "      <div class='partition p-cloud'><h3>Cloud (&lt;15%)</h3><div class='count'>{}</div></div>", partitions.cloud)?;
        writeln!(writer, "    </div>")?;
        writeln!(writer, "  </div>")?;

        // U-shape gene frequency bar chart
        writeln!(writer, "  <div class='section'>")?;
        writeln!(writer, "    <h2>Gene Frequency Distribution</h2>")?;
        writeln!(writer, "    <div class='chart' id='freq-chart'></div>")?;
        writeln!(writer, "  </div>")?;

        // Rarefaction curve
        writeln!(writer, "  <div class='section'>")?;
        writeln!(writer, "    <h2>Pangenome Rarefaction Curve</h2>")?;
        writeln!(writer, "    <div class='chart' id='rarefaction-chart'></div>")?;
        writeln!(writer, "  </div>")?;

        // D3.js rendering
        writeln!(writer, "  <script>")?;

        // Gene frequency bar chart
        writeln!(writer, "    const freqData = {};", freq_data_json)?;
        writeln!(writer, "    const freqChart = document.getElementById('freq-chart');")?;
        writeln!(writer, "    const fw = freqChart.clientWidth || 700;")?;
        writeln!(writer, "    const fh = 300;")?;
        writeln!(writer, "    const freqSvg = d3.select(freqChart).append('svg').attr('width', fw).attr('height', fh);")?;
        writeln!(writer, "    const xFreq = d3.scaleBand().domain(freqData.map(d => d[0])).range([60, fw-20]).padding(0.1);")?;
        writeln!(writer, "    const yMax = d3.max(freqData, d => d[1]) || 1;")?;
        writeln!(writer, "    const yFreq = d3.scaleLinear().domain([0, yMax * 1.1]).range([fh-30, 10]);")?;
        writeln!(writer, "    freqSvg.selectAll('.bar-freq').data(freqData).enter().append('rect')")?;
        writeln!(writer, "      .attr('class', 'bar-freq').attr('x', d => xFreq(d[0])).attr('y', d => yFreq(d[1]))")?;
        writeln!(writer, "      .attr('width', xFreq.bandwidth()).attr('height', d => fh - 30 - yFreq(d[1]));")?;
        writeln!(writer, "    freqSvg.append('g').attr('transform', 'translate(0,' + (fh-30) + ')').call(d3.axisBottom(xFreq).tickValues(xFreq.domain().filter((d,i) => !(i % Math.max(1, Math.floor(xFreq.domain().length/20))))));")?;
        writeln!(writer, "    freqSvg.append('g').attr('transform', 'translate(60,0)').call(d3.axisLeft(yFreq).ticks(6));")?;
        writeln!(writer, "    freqSvg.append('text').attr('class', 'axis-label').attr('x', fw/2).attr('y', fh-2).attr('text-anchor', 'middle').text('# Genomes');")?;
        writeln!(writer, "    freqSvg.append('text').attr('class', 'axis-label').attr('transform', 'rotate(-90)').attr('x', -fh/2).attr('y', 18).attr('text-anchor', 'middle').text('# Gene Families');")?;

        // Rarefaction line chart
        writeln!(writer, "    const rarData = {};", rarefaction_json)?;
        writeln!(writer, "    const rarChart = document.getElementById('rarefaction-chart');")?;
        writeln!(writer, "    const rw = rarChart.clientWidth || 700;")?;
        writeln!(writer, "    const rh = 300;")?;
        writeln!(writer, "    const rarSvg = d3.select(rarChart).append('svg').attr('width', rw).attr('height', rh);")?;
        writeln!(writer, "    const xRar = d3.scaleLinear().domain([0, d3.max(rarData, d => d[0]) || 1]).range([60, rw-20]);")?;
        writeln!(writer, "    const yRar = d3.scaleLinear().domain([0, d3.max(rarData, d => d[1]) || 1]).range([rh-30, 10]);")?;
        writeln!(writer, "    const rarLine = d3.line().x(d => xRar(d[0])).y(d => yRar(d[1]));")?;
        writeln!(writer, "    rarSvg.append('path').datum(rarData).attr('class', 'rarefaction-line').attr('d', rarLine);")?;
        writeln!(writer, "    rarSvg.selectAll('.rarefaction-point').data(rarData).enter().append('circle')")?;
        writeln!(writer, "      .attr('class', 'rarefaction-point').attr('cx', d => xRar(d[0])).attr('cy', d => yRar(d[1])).attr('r', 4);")?;
        writeln!(writer, "    rarSvg.append('g').attr('transform', 'translate(0,' + (rh-30) + ')').call(d3.axisBottom(xRar).ticks(8));")?;
        writeln!(writer, "    rarSvg.append('g').attr('transform', 'translate(60,0)').call(d3.axisLeft(yRar).ticks(6));")?;
        writeln!(writer, "    rarSvg.append('text').attr('class', 'axis-label').attr('x', rw/2).attr('y', rh-2).attr('text-anchor', 'middle').text('# Genomes Added');")?;
        writeln!(writer, "    rarSvg.append('text').attr('class', 'axis-label').attr('transform', 'rotate(-90)').attr('x', -rh/2).attr('y', 18).attr('text-anchor', 'middle').text('Cumulative Pangenome Size');")?;

        writeln!(writer, "  </script>")?;
        writeln!(writer, "</body>")?;
        writeln!(writer, "</html>")?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_partitions_basic() {
        // 10 genomes total
        let counts = vec![
            (10, 50),  // 10/10 = 100% -> core
            (9, 30),   // 9/10 = 90% -> shell
            (8, 20),   // 8/10 = 80% -> shell
            (5, 40),   // 5/10 = 50% -> shell
            (1, 100),  // 1/10 = 10% -> cloud
        ];
        let partitions = AbundanceVizWriter::compute_partitions(&counts, 10);
        assert_eq!(partitions.core, 50);
        assert_eq!(partitions.soft_core, 0);
        assert_eq!(partitions.shell, 90); // 30 + 20 + 40
        assert_eq!(partitions.cloud, 100);
    }

    #[test]
    fn test_compute_partitions_soft_core() {
        // 100 genomes total
        let counts = vec![
            (100, 200), // 100% -> core
            (98, 30),   // 98% -> soft-core
            (96, 25),   // 96% -> soft-core
            (50, 60),   // 50% -> shell
            (10, 80),   // 10% -> cloud
        ];
        let partitions = AbundanceVizWriter::compute_partitions(&counts, 100);
        assert_eq!(partitions.core, 200);
        assert_eq!(partitions.soft_core, 55); // 30 + 25
        assert_eq!(partitions.shell, 60);
        assert_eq!(partitions.cloud, 80);
    }

    #[test]
    fn test_compute_partitions_empty() {
        let partitions = AbundanceVizWriter::compute_partitions(&[], 10);
        assert_eq!(partitions.core, 0);
        assert_eq!(partitions.soft_core, 0);
        assert_eq!(partitions.shell, 0);
        assert_eq!(partitions.cloud, 0);
    }

    #[test]
    fn test_compute_partitions_zero_genomes() {
        let counts = vec![(0, 10), (1, 20)];
        let partitions = AbundanceVizWriter::compute_partitions(&counts, 0);
        assert_eq!(partitions.core, 0);
        assert_eq!(partitions.soft_core, 0);
        assert_eq!(partitions.shell, 0);
        assert_eq!(partitions.cloud, 0);
    }

    #[test]
    fn test_compute_partitions_boundary_95_percent() {
        // 20 genomes: 19/20 = 95% exactly -> soft-core
        let counts = vec![(19, 10)];
        let partitions = AbundanceVizWriter::compute_partitions(&counts, 20);
        assert_eq!(partitions.soft_core, 10);
        assert_eq!(partitions.shell, 0);
    }

    #[test]
    fn test_compute_partitions_boundary_15_percent() {
        // 20 genomes: 3/20 = 15% exactly -> shell
        let counts = vec![(3, 10)];
        let partitions = AbundanceVizWriter::compute_partitions(&counts, 20);
        assert_eq!(partitions.shell, 10);
        assert_eq!(partitions.cloud, 0);
    }

    #[test]
    fn test_write_report_creates_html() {
        let presence_counts = vec![
            (10, 50),
            (5, 30),
            (1, 100),
        ];
        let rarefaction_data = vec![
            (1, 500),
            (5, 1200),
            (10, 1800),
        ];
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("abundance_report.html");

        AbundanceVizWriter::write_report(
            &presence_counts,
            &rarefaction_data,
            10,
            180,
            &path,
        ).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();

        // Key structural elements
        assert!(content.contains("PanMiner Gene Abundance Report"));
        assert!(content.contains("d3.v7.min.js"));
        assert!(content.contains("freq-chart"));
        assert!(content.contains("rarefaction-chart"));

        // Data embedded as JSON
        assert!(content.contains("freqData"));
        assert!(content.contains("rarData"));

        // Partition values
        assert!(content.contains("50</div>"));  // core count
        assert!(content.contains("100</div>")); // cloud count

        // SVG elements are created by D3.js at runtime, but the chart divs exist
        assert!(content.contains("id='freq-chart'"));
        assert!(content.contains("id='rarefaction-chart'"));
    }

    #[test]
    fn test_write_report_contains_partition_labels() {
        let presence_counts = vec![(10, 1)];
        let rarefaction_data = vec![(1, 1)];
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("abundance_report.html");

        AbundanceVizWriter::write_report(
            &presence_counts,
            &rarefaction_data,
            10,
            1,
            &path,
        ).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("Core"));
        assert!(content.contains("Soft-core"));
        assert!(content.contains("Shell"));
        assert!(content.contains("Cloud"));
    }

    #[test]
    fn test_write_report_empty_data() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("abundance_report.html");

        AbundanceVizWriter::write_report(
            &[],
            &[],
            0,
            0,
            &path,
        ).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("PanMiner Gene Abundance Report"));
        assert!(content.contains("d3.v7.min.js"));
    }
}