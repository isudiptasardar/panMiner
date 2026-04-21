//! Pre-QC diagnostic HTML report module.
//!
//! Generates an interactive HTML report with MDS scatter plots and
//! contamination/completeness bar charts from QC results.

use crate::error::Result;
use crate::io::{GenomeQC, MdsProjection};
use std::path::Path;
use std::io::Write;

/// Escape special HTML characters to prevent XSS injection.
fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#x27;")
}

/// Write a pre-QC HTML report with MDS scatter plot and contamination bar chart.
pub fn write_qc_html_report(
    qc_results: &[GenomeQC],
    mds: Option<&MdsProjection>,
    path: &Path,
) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);

    writeln!(writer, "<!DOCTYPE html>")?;
    writeln!(writer, "<html lang='en'>")?;
    writeln!(writer, "<head>")?;
    writeln!(writer, "  <meta charset='UTF-8'>")?;
    writeln!(writer, "  <meta name='viewport' content='width=device-width, initial-scale=1.0'>")?;
    writeln!(writer, "  <title>PanMiner QC Report</title>")?;
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
    writeln!(writer, "    .bar-pass {{ fill: #4ecdc4; }}")?;
    writeln!(writer, "    .bar-fail {{ fill: #e94560; }}")?;
    writeln!(writer, "    .point {{ fill: #e94560; cursor: pointer; opacity: 0.8; }}")?;
    writeln!(writer, "    .point:hover {{ opacity: 1; fill: #4ecdc4; }}")?;
    writeln!(writer, "    .axis-label {{ fill: #666; font-size: 12px; }}")?;
    writeln!(writer, "    .tooltip {{ position: absolute; background: rgba(0,0,0,0.95); padding: 10px 14px; border-radius: 6px; font-size: 13px; pointer-events: none; border: 1px solid rgba(255,255,255,0.1); }}")?;
    writeln!(writer, "    .tooltip h4 {{ color: #e94560; margin-bottom: 6px; }}")?;
    writeln!(writer, "    .tooltip .meta {{ color: #aaa; font-size: 12px; line-height: 1.6; }}")?;
    writeln!(writer, "    table {{ width: 100%; border-collapse: collapse; }}")?;
    writeln!(writer, "    th, td {{ padding: 10px 12px; text-align: left; border-bottom: 1px solid rgba(255,255,255,0.05); }}")?;
    writeln!(writer, "    th {{ color: #888; font-size: 12px; text-transform: uppercase; letter-spacing: 0.5px; }}")?;
    writeln!(writer, "    .pass {{ color: #4ecdc4; }}")?;
    writeln!(writer, "    .fail {{ color: #e94560; }}")?;
    writeln!(writer, "  </style>")?;
    writeln!(writer, "</head>")?;
    writeln!(writer, "<body>")?;

    // Header with stats
    let total = qc_results.len();
    let passed: usize = qc_results.iter().filter(|q| q.passed).count();
    let failed = total - passed;
    let avg_comp = if total > 0 { qc_results.iter().map(|q| q.completeness).sum::<f64>() / total as f64 } else { 0.0 };
    let avg_cont = if total > 0 { qc_results.iter().map(|q| q.contamination).sum::<f64>() / total as f64 } else { 0.0 };

    writeln!(writer, "  <div class='header'>")?;
    writeln!(writer, "    <h1>PanMiner QC Report</h1>")?;
    writeln!(writer, "    <div class='stats'>")?;
    writeln!(writer, "      <div class='stat-item'><div class='stat-value'>{}</div><div class='stat-label'>Total Genomes</div></div>", total)?;
    writeln!(writer, "      <div class='stat-item'><div class='stat-value'>{}</div><div class='stat-label'>Passed</div></div>", passed)?;
    writeln!(writer, "      <div class='stat-item'><div class='stat-value'>{}</div><div class='stat-label'>Failed</div></div>", failed)?;
    writeln!(writer, "      <div class='stat-item'><div class='stat-value'>{:.1}%</div><div class='stat-label'>Avg Completeness</div></div>", avg_comp)?;
    writeln!(writer, "      <div class='stat-item'><div class='stat-value'>{:.1}%</div><div class='stat-label'>Avg Contamination</div></div>", avg_cont)?;
    writeln!(writer, "    </div>")?;
    writeln!(writer, "  </div>")?;

    // MDS Scatter Plot
    if let Some(_mds_data) = mds {
        writeln!(writer, "  <div class='section'>")?;
        writeln!(writer, "    <h2>Genomic Distance MDS Projection</h2>")?;
        writeln!(writer, "    <div class='chart' id='mds-chart'></div>")?;
        writeln!(writer, "  </div>")?;
    }

    // Contamination Bar Chart
    writeln!(writer, "  <div class='section'>")?;
    writeln!(writer, "    <h2>Assembly Completeness</h2>")?;
    writeln!(writer, "    <div class='chart' id='comp-chart'></div>")?;
    writeln!(writer, "  </div>")?;

    writeln!(writer, "  <div class='section'>")?;
    writeln!(writer, "    <h2>Assembly Contamination</h2>")?;
    writeln!(writer, "    <div class='chart' id='cont-chart'></div>")?;
    writeln!(writer, "  </div>")?;

    // Data table
    writeln!(writer, "  <div class='section'>")?;
    writeln!(writer, "    <h2>Genome Details</h2>")?;
    writeln!(writer, "    <table>")?;
    writeln!(writer, "      <thead><tr><th>Genome</th><th>Completeness</th><th>Contamination</th><th>Status</th></tr></thead>")?;
    writeln!(writer, "      <tbody>")?;
    for qc in qc_results {
        let status = if qc.passed { "PASS" } else { "FAIL" };
        let status_class = if qc.passed { "pass" } else { "fail" };
        writeln!(writer, "        <tr>")?;
        writeln!(writer, "          <td>{}</td>", html_escape(&qc.genome_id))?;
        writeln!(writer, "          <td>{:.1}%</td>", qc.completeness)?;
        writeln!(writer, "          <td>{:.1}%</td>", qc.contamination)?;
        writeln!(writer, "          <td class='{}'>{}</td>", status_class, status)?;
        writeln!(writer, "        </tr>")?;
    }
    writeln!(writer, "      </tbody>")?;
    writeln!(writer, "    </table>")?;
    writeln!(writer, "  </div>")?;

    writeln!(writer, "  <div class='tooltip' id='tooltip' style='display:none;'></div>")?;

    // D3.js charts
    writeln!(writer, "  <script>")?;
    writeln!(writer, "    function escapeHtml(s) {{ return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/\"/g,'&quot;').replace(/'/g,'&#x27;'); }}")?;
    writeln!(writer, "    const tooltip = document.getElementById('tooltip');")?;

    // Completeness bar chart
    writeln!(writer, "    const qcData = {};", serde_json::to_string(qc_results).unwrap_or_else(|_| "[]".to_string()))?;
    writeln!(writer, "    const compChart = document.getElementById('comp-chart');")?;
    writeln!(writer, "    const w = compChart.clientWidth || 700;")?;
    writeln!(writer, "    const h = 200;")?;
    writeln!(writer, "    const compSvg = d3.select(compChart).append('svg').attr('width', w).attr('height', h);")?;
    writeln!(writer, "    const xComp = d3.scaleBand().domain(qcData.map((_,i) => i)).range([50, w-10]).padding(0.2);")?;
    writeln!(writer, "    const yComp = d3.scaleLinear().domain([0, 100]).range([h-20, 10]);")?;
    writeln!(writer, "    compSvg.selectAll('.bar').data(qcData).enter().append('rect')")?;
    writeln!(writer, "      .attr('class', d => d.passed ? 'bar-pass' : 'bar-fail').attr('x', (_,i) => xComp(i))")?;
    writeln!(writer, "      .attr('y', d => yComp(d.completeness)).attr('width', xComp.bandwidth())")?;
    writeln!(writer, "      .attr('height', d => h - 20 - yComp(d.completeness))")?;
    writeln!(writer, "      .on('mouseover', (event, d) => {{ tooltip.style.display='block'; tooltip.innerHTML='<h4>'+escapeHtml(d.genome_id)+'</h4><div class=meta>Completeness: '+d.completeness.toFixed(1)+'%</div>'; tooltip.style.left=(event.pageX+10)+'px'; tooltip.style.top=(event.pageY-10)+'px'; }})")?;
    writeln!(writer, "      .on('mouseout', () => {{ tooltip.style.display='none'; }});")?;
    writeln!(writer, "    compSvg.append('text').attr('class', 'axis-label').attr('x', w/2).attr('y', h-2).attr('text-anchor', 'middle').text('Genome');")?;
    writeln!(writer, "    compSvg.append('text').attr('class', 'axis-label').attr('transform', 'rotate(-90)').attr('x', -h/2).attr('y', 15).attr('text-anchor', 'middle').text('Completeness (%)');")?;

    // Contamination bar chart
    writeln!(writer, "    const contChart = document.getElementById('cont-chart');")?;
    writeln!(writer, "    const contSvg = d3.select(contChart).append('svg').attr('width', w).attr('height', h);")?;
    writeln!(writer, "    const xCont = d3.scaleBand().domain(qcData.map((_,i) => i)).range([50, w-10]).padding(0.2);")?;
    writeln!(writer, "    const maxCont = Math.min(50, (d3.max(qcData, d => d.contamination) || 0) * 1.1 || 1);")?;
    writeln!(writer, "    const yCont = d3.scaleLinear().domain([0, maxCont]).range([h-20, 10]);")?;
    writeln!(writer, "    contSvg.selectAll('.bar').data(qcData).enter().append('rect')")?;
    writeln!(writer, "      .attr('class', d => d.passed ? 'bar-pass' : 'bar-fail').attr('x', (_,i) => xCont(i))")?;
    writeln!(writer, "      .attr('y', d => yCont(d.contamination)).attr('width', xCont.bandwidth())")?;
    writeln!(writer, "      .attr('height', d => h - 20 - yCont(d.contamination))")?;
    writeln!(writer, "      .on('mouseover', (event, d) => {{ tooltip.style.display='block'; tooltip.innerHTML='<h4>'+escapeHtml(d.genome_id)+'</h4><div class=meta>Contamination: '+d.contamination.toFixed(1)+'%</div>'; tooltip.style.left=(event.pageX+10)+'px'; tooltip.style.top=(event.pageY-10)+'px'; }})")?;
    writeln!(writer, "      .on('mouseout', () => {{ tooltip.style.display='none'; }});")?;
    writeln!(writer, "    contSvg.append('text').attr('class', 'axis-label').attr('x', w/2).attr('y', h-2).attr('text-anchor', 'middle').text('Genome');")?;
    writeln!(writer, "    contSvg.append('text').attr('class', 'axis-label').attr('transform', 'rotate(-90)').attr('x', -h/2).attr('y', 15).attr('text-anchor', 'middle').text('Contamination (%)');")?;

    // MDS scatter plot
    if let Some(mds_data) = mds {
        let coords_json = serde_json::to_string(&mds_data.coordinates).unwrap_or_else(|_| "[]".to_string());
        let labels_json = serde_json::to_string(&mds_data.labels).unwrap_or_else(|_| "[]".to_string());
        writeln!(writer, "    const mdsChart = document.getElementById('mds-chart');")?;
        writeln!(writer, "    const mdsW = mdsChart.clientWidth || 700;")?;
        writeln!(writer, "    const mdsH = 350;")?;
        writeln!(writer, "    const mdsSvg = d3.select(mdsChart).append('svg').attr('width', mdsW).attr('height', mdsH);")?;
        writeln!(writer, "    const mdsCoords = {};", coords_json)?;
        writeln!(writer, "    const mdsLabels = {};", labels_json)?;
        writeln!(writer, "    const xMds = d3.scaleLinear().domain(d3.extent(mdsCoords, d => d[0])).range([50, mdsW-20]);")?;
        writeln!(writer, "    const yMds = d3.scaleLinear().domain(d3.extent(mdsCoords, d => d[1])).range([mdsH-20, 10]);")?;
        writeln!(writer, "    mdsSvg.selectAll('.point').data(mdsCoords).enter().append('circle')")?;
        writeln!(writer, "      .attr('class', 'point').attr('cx', d => xMds(d[0])).attr('cy', d => yMds(d[1])).attr('r', 7)")?;
        writeln!(writer, "      .on('mouseover', (event, d) => {{ const i = mdsCoords.indexOf(d); tooltip.style.display='block'; tooltip.innerHTML='<h4>'+escapeHtml(mdsLabels[i])+'</h4>'; tooltip.style.left=(event.pageX+10)+'px'; tooltip.style.top=(event.pageY-10)+'px'; }})")?;
        writeln!(writer, "      .on('mouseout', () => {{ tooltip.style.display='none'; }});")?;
        writeln!(writer, "    mdsSvg.append('text').attr('class', 'axis-label').attr('x', mdsW/2).attr('y', mdsH-2).attr('text-anchor', 'middle').text('MDS Dimension 1');")?;
        writeln!(writer, "    mdsSvg.append('text').attr('class', 'axis-label').attr('transform', 'rotate(-90)').attr('x', -mdsH/2).attr('y', 15).attr('text-anchor', 'middle').text('MDS Dimension 2');")?;
    }

    writeln!(writer, "  </script>")?;
    writeln!(writer, "</body>")?;
    writeln!(writer, "</html>")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_qc_html_report() {
        let qc_results = vec![
            GenomeQC { genome_id: "g1".to_string(), completeness: 95.0, contamination: 2.0, passed: true, ..Default::default() },
            GenomeQC { genome_id: "g2".to_string(), completeness: 60.0, contamination: 15.0, passed: false, ..Default::default() },
        ];
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("qc_report.html");
        write_qc_html_report(&qc_results, None, &path).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("PanMiner QC Report"));
        assert!(content.contains("Total Genomes"));
        assert!(content.contains("g1"));
        assert!(content.contains("g2"));
    }

    #[test]
    fn test_write_qc_html_report_with_mds() {
        let qc_results = vec![
            GenomeQC { genome_id: "g1".to_string(), completeness: 95.0, contamination: 2.0, passed: true, ..Default::default() },
        ];
        let mds = MdsProjection {
            coordinates: vec![(0.1, 0.2)],
            labels: vec!["g1".to_string()],
        };
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("qc_report.html");
        write_qc_html_report(&qc_results, Some(&mds), &path).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("mds-chart"));
    }
}
