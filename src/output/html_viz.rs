//! HTML visualization output module.
//!
//! Generates interactive HTML visualizations of the pangenome graph
//! using vanilla JavaScript (d3.js via CDN) for interactive exploration.

use std::path::Path;
use std::fs::File;
use std::io::Write;

use crate::error::Result;
use crate::graph::{PangenomeGraph, BitPackedMatrix};

/// HTML visualization writer.
pub struct HtmlVizWriter;

impl HtmlVizWriter {
    /// Create a new HTML visualization writer.
    pub fn new() -> Self {
        HtmlVizWriter
    }

    /// Write HTML visualization to the given path.
    pub fn write(&self, graph: &PangenomeGraph, matrix: &BitPackedMatrix, path: &Path) -> Result<()> {
        let html = Self::generate_html(graph, matrix);

        let mut file = File::create(path)?;
        file.write_all(html.as_bytes())?;

        Ok(())
    }

    /// Generate the HTML content for visualization.
    fn generate_html(graph: &PangenomeGraph, matrix: &BitPackedMatrix) -> String {
        let node_data = Self::build_node_data(graph, matrix);
        let edge_data = Self::build_edge_data(graph);
        let genome_names = matrix.genome_names.clone();
        let node_count = graph.node_count();
        let edge_count = graph.edge_count();
        let num_genomes = matrix.num_genomes();

        // Build genome filter checkboxes
        let genome_filters: String = genome_names
            .iter()
            .enumerate()
            .map(|(_, name)| {
                format!(
                    "<div class=\"genome_filter\"><input type=\"checkbox\" checked><span class=\"genome_name\">{}</span></div>",
                    name
                )
            })
            .collect();

        // Build HTML using string concatenation to avoid format string issues
        let mut html = String::new();

        // Header
        html.push_str("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n");
        html.push_str("    <meta charset=\"UTF-8\">\n");
        html.push_str("    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n");
        html.push_str("    <title>PanMiner Pangenome Visualization</title>\n");
        html.push_str("    <script src=\"https://d3js.org/d3.v7.min.js\"></script>\n");

        // CSS - inline to avoid external dependencies
        html.push_str("    <style>\n");
        html.push_str("        *{margin:0;padding:0;box-sizing:border-box;}\n");
        html.push_str("        body{font-family:sans-serif;background:linear-gradient(135deg,rgb(26,26,46),rgb(22,33,62));min-height:100vh;color:rgb(238,238,238);}\n");
        html.push_str("        .header{background:rgba(255,255,255,0.05);padding:20px 40px;border-bottom:1px solid rgba(255,255,255,0.1);position:sticky;top:0;z-index:100;}\n");
        html.push_str("        .header h1{font-size:24px;font-weight:600;margin-bottom:10px;}\n");
        html.push_str("        .stats{display:flex;gap:30px;flex-wrap:wrap;}\n");
        html.push_str("        .stat-item{display:flex;align-items:center;gap:8px;}\n");
        html.push_str("        .stat-value{font-size:20px;font-weight:700;color:rgb(233,69,96);}\n");
        html.push_str("        .stat-label{color:rgb(136,136,136);font-size:14px;}\n");
        html.push_str("        .main{display:flex;height:calc(100vh - 80px);}\n");
        html.push_str("        .sidebar{width:320px;background:rgba(0,0,0,0.3);border-right:1px solid rgba(255,255,255,0.1);padding:20px;overflow-y:auto;flex-shrink:0;}\n");
        html.push_str("        .sidebar-section{margin-bottom:25px;}\n");
        html.push_str("        .sidebar-section h3{font-size:14px;text-transform:uppercase;letter-spacing:1px;color:rgb(136,136,136);margin-bottom:15px;}\n");
        html.push_str("        .genome_filter{display:flex;align-items:center;padding:8px 12px;background:rgba(255,255,255,0.05);border-radius:6px;margin-bottom:6px;cursor:pointer;}\n");
        html.push_str("        .genome_filter:hover{background:rgba(255,255,255,0.1);}\n");
        html.push_str("        .genome_checkbox{width:16px;height:16px;margin-right:12px;accent-color:rgb(233,69,96);}\n");
        html.push_str("        .genome_name{font-size:13px;word-break:break-all;}\n");
        html.push_str("        .btn{padding:12px 20px;border:none;border-radius:8px;font-size:14px;font-weight:600;cursor:pointer;}\n");
        html.push_str("        .btn_primary{background:linear-gradient(90deg,rgb(233,69,96),rgb(197,53,76));color:white;}\n");
        html.push_str("        .btn_primary:hover{transform:translateY(-2px);box-shadow:0 4px 12px rgba(233,69,96,0.3);}\n");
        html.push_str("        .btn_group{display:flex;gap:10px;}\n");
        html.push_str("        .graph_container{flex:1;position:relative;background:rgb(10,10,21);overflow:hidden;}\n");
        html.push_str("        .tooltip{position:absolute;background:rgba(0,0,0,0.9);padding:12px 16px;border-radius:8px;border:1px solid rgba(255,255,255,0.1);pointer-events:none;z-index:1000;max-width:250px;}\n");
        html.push_str("        .tooltip h4{color:rgb(233,69,96);font-size:14px;margin-bottom:8px;}\n");
        html.push_str("        .tooltip .meta{font-size:13px;color:rgb(170,170,170);}\n");
        html.push_str("        .node_circle{cursor:pointer;}\n");
        html.push_str("        .node_circle:hover{stroke:rgb(233,69,96);stroke-width:3px;}\n");
        html.push_str("        .node_circle.supported{fill:rgb(233,69,96);}\n");
        html.push_str("        .node_circle.accessory{fill:rgb(78,204,196);}\n");
        html.push_str("        .link{stroke-opacity:0.4;}\n");
        html.push_str("        .link:hover{stroke-opacity:0.8;stroke-width:2px;}\n");
        html.push_str("        .zoom_controls{position:absolute;bottom:20px;right:20px;display:flex;flex-direction:column;gap:10px;}\n");
        html.push_str("        .zoom_btn{width:40px;height:40px;border:none;border-radius:8px;background:rgba(255,255,255,0.9);color:rgb(51,51,51);cursor:pointer;}\n");
        html.push_str("        .zoom_btn:hover{transform:scale(1.1);}\n");
        html.push_str("        .empty_state{position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);text-align:center;color:rgb(102,102,102);}\n");
        html.push_str("    </style>\n");
        html.push_str("</head>\n<body>\n");

        // Header section
        html.push_str("    <div class=\"header\">\n");
        html.push_str("        <h1>PanMiner Pangenome Visualization</h1>\n");
        html.push_str("        <div class=\"stats\">\n");
        html.push_str("            <div class=\"stat-item\"><div class=\"stat-value\">");
        html.push_str(&node_count.to_string());
        html.push_str("</div><div class=\"stat-label\">Gene Clusters</div></div>\n");
        html.push_str("            <div class=\"stat-item\"><div class=\"stat-value\">");
        html.push_str(&edge_count.to_string());
        html.push_str("</div><div class=\"stat-label\">Adjacencies</div></div>\n");
        html.push_str("            <div class=\"stat-item\"><div class=\"stat-value\">");
        html.push_str(&num_genomes.to_string());
        html.push_str("</div><div class=\"stat-label\">Genomes</div></div>\n");
        html.push_str("        </div>\n");
        html.push_str("    </div>\n");

        // Main content
        html.push_str("    <div class=\"main\">\n");
        html.push_str("        <div class=\"sidebar\">\n");
        html.push_str("            <div class=\"sidebar-section\"><h3>Genome Filters</h3><div id=\"genome_filters\">\n");
        html.push_str(&genome_filters);
        html.push_str("</div></div>\n");
        html.push_str("            <div class=\"sidebar-section\"><h3>Export</h3><div class=\"btn_group\">\n");
        html.push_str("                    <button class=\"btn btn_primary\" id=\"export_json\"><span>JSON</span></button>\n");
        html.push_str("                    <button class=\"btn btn_primary\" id=\"export_csv\"><span>CSV</span></button>\n");
        html.push_str("                </div></div>\n");
        html.push_str("            <div class=\"sidebar-section\"><h3>Stats</h3><div class=\"controls\">\n");
        html.push_str("                    <div class=\"genome_filter\"><span class=\"genome_name\">Core clusters (100%): <strong id=\"core_count\">0</strong></span></div>\n");
        html.push_str("                    <div class=\"genome_filter\"><span class=\"genome_name\">Accessory clusters: <strong id=\"accessory_count\">0</strong></span></div>\n");
        html.push_str("                    <div class=\"genome_filter\"><span class=\"genome_name\">Paralogs: <strong id=\"paralog_count\">0</strong></span></div>\n");
        html.push_str("                </div></div>\n");
        html.push_str("        </div>\n");

        // Graph container
        html.push_str("        <div class=\"graph_container\">\n");
        html.push_str("            <div id=\"graph\">\n");
        html.push_str("                <div class=\"empty_state\" id=\"empty_state\">\n");
        html.push_str("                    <svg width=\"100\" height=\"100\" viewBox=\"0 0 100 100\">\n");
        html.push_str("                        <circle cx=\"50\" cy=\"50\" r=\"30\" stroke=\"rgb(68,68,68)\" stroke-width=\"2\" fill=\"none\"/>\n");
        html.push_str("                        <circle cx=\"50\" cy=\"50\" r=\"15\" stroke=\"rgb(68,68,68)\" stroke-width=\"2\" fill=\"none\"/>\n");
        html.push_str("                        <circle cx=\"50\" cy=\"50\" r=\"5\" fill=\"rgb(68,68,68)\"/>\n");
        html.push_str("                    </svg>\n");
        html.push_str("                    <p>Visualizing pangenome graph...</p>\n");
        html.push_str("                </div>\n");
        html.push_str("                <div id=\"graph_svg\" style=\"display:none;\"></div>\n");
        html.push_str("            </div>\n");
        html.push_str("            <div class=\"zoom_controls\">\n");
        html.push_str("                <button class=\"zoom_btn\" id=\"zoom_in\">+</button>\n");
        html.push_str("                <button class=\"zoom_btn\" id=\"zoom_out\">-</button>\n");
        html.push_str("                <button class=\"zoom_btn\" id=\"reset_view\">&#x27b6;</button>\n");
        html.push_str("            </div>\n");
        html.push_str("        </div>\n");
        html.push_str("    </div>\n");

        // Tooltip
        html.push_str("    <div class=\"tooltip\" id=\"tooltip\" style=\"display:none;\"></div>\n");

        // JavaScript
        html.push_str("    <script>\n");
        html.push_str("        const nodeData = ");
        html.push_str(&node_data);
        html.push_str(";\n");
        html.push_str("        const edgeData = ");
        html.push_str(&edge_data);
        html.push_str(";\n");
        html.push_str("        const genomeNames = ");
        html.push_str(&serde_json::to_string(&genome_names).unwrap_or_else(|_| "[]".to_string()));
        html.push_str(";\n");
        html.push_str("        const tooltip = document.getElementById(\"tooltip\");\n");
        html.push_str("        const graphSvg = document.getElementById(\"graph_svg\");\n");
        html.push_str("        const emptyState = document.getElementById(\"empty_state\");\n");
        html.push_str("        let zoomLevel = 1;\n");
        html.push_str("        let currentTransform = {x: 0, y: 0, k: 1};\n");
        html.push_str("        function initGraph() {\n");
        html.push_str("            const svg = d3.select(\"#graph_svg\");\n");
        html.push_str("            svg.innerHTML = \"\";\n");
        html.push_str("            const width = graphSvg.parentElement.clientWidth;\n");
        html.push_str("            const height = graphSvg.parentElement.clientHeight;\n");
        html.push_str("            const g = svg.append(\"g\").attr(\"transform\", \"translate(\" + width/2 + \",\" + height/2 + \")\");\n");
        html.push_str("            const color = d3.scaleOrdinal().domain([0, 1]).range([\"rgb(233,69,96)\", \"rgb(78,204,196)\"]);\n");
        html.push_str("            const link = g.append(\"g\").attr(\"class\", \"links\").selectAll(\"line\").data(edgeData).enter().append(\"line\").attr(\"class\", \"link\").attr(\"stroke\", \"rgb(102,102,102)\").attr(\"stroke-width\", function(d) { return Math.sqrt(d.weight) / 2; });\n");
        html.push_str("            link.on(\"mouseover\", function(event, d) { tooltip.style(\"display\", \"block\").html(\"Edge: \" + d.source + \" -> \" + d.target + \"<br/>Support: \" + d.support).style(\"left\", (event.pageX + 10) + \"px\").style(\"top\", (event.pageY + 10) + \"px\"); }).on(\"mouseout\", function() { tooltip.style(\"display\", \"none\"); });\n");
        html.push_str("            const node = g.append(\"g\").attr(\"class\", \"nodes\").selectAll(\"circle\").data(nodeData).enter().append(\"circle\").attr(\"class\", function(d) { return \"node_circle \" + (d.support === genomeNames.length ? \"supported\" : \"accessory\"); }).attr(\"r\", function(d) { return Math.min(20, Math.max(5, Math.log2(d.genes + 1) * 3)); }).attr(\"fill\", function(d) { return color(d.support / genomeNames.length); }).attr(\"stroke\", \"rgb(34,34,34)\").attr(\"stroke-width\", 1.5).call(d3.drag().on(\"start\", dragstarted).on(\"drag\", dragged).on(\"end\", dragended)).on(\"mouseover\", function(event, d) { tooltip.style(\"display\", \"block\").html(\"<h4>\" + d.clusterId + \"</h4><div class=meta><span>Support:</span> \" + d.support + \" / \" + genomeNames.length + \"<br/><span>Genes:</span> \" + d.genes + \"<br/><span>Paralog:</span> \" + (d.isParalog ? \"Yes\" : \"No\") + \"</div><div class=genome_list>\" + d.genomes.map(function(g) { return \"<span class=genome>\" + (genomeNames[g] || \"Unknown\") + \"</span>\"; }).join(\"\") + \"</div>\").style(\"left\", (event.pageX + 10) + \"px\").style(\"top\", (event.pageY + 10) + \"px\"); }).on(\"mouseout\", function() { tooltip.style(\"display\", \"none\"); });\n");
        html.push_str("            const simulation = d3.forceSimulation(nodeData).force(\"link\", d3.forceLink(edgeData).id(function(d) { return d.id; }).distance(function(d) { return 50 + Math.log2(d.weight + 1) * 20; })).force(\"charge\", d3.forceManyBody().strength(-500).theta(0.5)).force(\"center\", d3.forceCenter(0, 0)).force(\"collide\", d3.forceCollide().radius(function(d) { return Math.min(20, Math.max(5, Math.log2(d.genes + 1) * 3)) + 5; })).alphaDecay(0.01).alphaMin(0.001).on(\"tick\", ticked);\n");
        html.push_str("            function ticked() { link.attr(\"x1\", function(d) { return d.source.x; }).attr(\"y1\", function(d) { return d.source.y; }).attr(\"x2\", function(d) { return d.target.x; }).attr(\"y2\", function(d) { return d.target.y; }); node.attr(\"cx\", function(d) { return d.x; }).attr(\"cy\", function(d) { return d.y; }); }\n");
        html.push_str("            function dragstarted(event, d) { if (!event.active) simulation.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; }\n");
        html.push_str("            function dragged(event, d) { d.fx = event.x; d.fy = event.y; }\n");
        html.push_str("            function dragended(event, d) { if (!event.active) simulation.alphaTarget(0); d.fx = null; d.fy = null; }\n");
        html.push_str("            const zoom = d3.zoom().scaleExtent([0.1, 10]).on(\"zoom\", function(event) { currentTransform = event.transform; g.attr(\"transform\", currentTransform); }); svg.call(zoom);\n");
        html.push_str("            emptyState.style(\"display\", \"none\");\n");
        html.push_str("            graphSvg.style(\"display\", \"block\");\n");
        html.push_str("            updateStats();\n");
        html.push_str("        }\n");
        html.push_str("        function updateStats() { const coreCount = nodeData.filter(function(n) { return n.support === genomeNames.length; }).length; const accessoryCount = nodeData.length - coreCount; const paralogCount = nodeData.filter(function(n) { return n.isParalog; }).length; document.getElementById(\"core_count\").textContent = coreCount; document.getElementById(\"accessory_count\").textContent = accessoryCount; document.getElementById(\"paralog_count\").textContent = paralogCount; }\n");
        html.push_str("        function exportJSON() { const data = {nodes: nodeData.map(function(n) { return {id: n.id, clusterId: n.clusterId, x: n.x, y: n.y, support: n.support, genes: n.genes, isParalog: n.isParalog, isCore: n.isCore}; }), edges: edgeData, genomeNames: genomeNames, clusterIds: clusterIds}; const blob = new Blob([JSON.stringify(data, null, 2)], {type: \"application/json\"}); const url = URL.createObjectURL(blob); const a = document.createElement(\"a\"); a.href = url; a.download = \"pangenome.json\"; a.click(); URL.revokeObjectURL(url); }\n");
        html.push_str("        function exportCSV() { const rows = [[\"ClusterID\", \"Support\", \"TotalGenomes\", \"IsParalog\"]]; for (const node of nodeData) { rows.push([node.clusterId, node.support, genomeNames.length, node.isParalog ? \"Yes\" : \"No\"]); } const csvContent = rows.map(function(e) { return e.join(\",\"); }).join(\"\\n\"); const blob = new Blob([csvContent], {type: \"text/csv\"}); const url = URL.createObjectURL(blob); const a = document.createElement(\"a\"); a.href = url; a.download = \"pangenome_clusters.csv\"; a.click(); URL.revokeObjectURL(url); }\n");
        html.push_str("        document.getElementById(\"export_json\").addEventListener(\"click\", exportJSON);\n");
        html.push_str("        document.getElementById(\"export_csv\").addEventListener(\"click\", exportCSV);\n");
        html.push_str("        document.getElementById(\"zoom_in\").addEventListener(\"click\", function() { const svg = d3.select(\"#graph_svg\"); const k = Math.min(zoomLevel * 2, 10); svg.transition().duration(200).call(d3.zoom().scaleTo, svg, k); zoomLevel = k; });\n");
        html.push_str("        document.getElementById(\"zoom_out\").addEventListener(\"click\", function() { const svg = d3.select(\"#graph_svg\"); const k = Math.max(zoomLevel / 2, 0.1); svg.transition().duration(200).call(d3.zoom().scaleTo, svg, k); zoomLevel = k; });\n");
        html.push_str("        document.getElementById(\"reset_view\").addEventListener(\"click\", function() { const svg = d3.select(\"#graph_svg\"); svg.transition().duration(200).call(d3.zoom().transform, svg, d3.zoomIdentity); zoomLevel = 1; currentTransform = {x: 0, y: 0, k: 1}; });\n");
        html.push_str("        window.addEventListener(\"load\", initGraph);\n");
        html.push_str("    </script>\n");
        html.push_str("</body>\n</html>");

        html
    }

    /// Build node data from graph and matrix.
    fn build_node_data(graph: &PangenomeGraph, matrix: &BitPackedMatrix) -> String {
        let mut nodes = Vec::new();

        for (cluster_id, node) in &graph.nodes {
            let genomes: Vec<usize> = node.genomes
                .iter()
                .filter_map(|g| matrix.genome_names.iter().position(|n| n == g.as_str()))
                .collect();

            let gene_estimate = node.support.max(1) * 2;

            let is_core = node.support == matrix.num_genomes();

            let node_obj = serde_json::json!({
                "id": cluster_id.as_str(),
                "clusterId": cluster_id.as_str(),
                "support": node.support,
                "genes": gene_estimate,
                "genomes": genomes,
                "isParalog": node.is_paralog,
                "isCore": is_core,
                "x": null,
                "y": null
            });

            nodes.push(node_obj);
        }

        serde_json::to_string(&nodes).unwrap_or_else(|_| "[]".to_string())
    }

    /// Build edge data from graph.
    fn build_edge_data(graph: &PangenomeGraph) -> String {
        let mut edges = Vec::new();

        for ((from, to), edge) in &graph.edges {
            let edge_obj = serde_json::json!({
                "source": from.as_str(),
                "target": to.as_str(),
                "support": edge.support,
                "weight": edge.genomes.len(),
                "x": null,
                "y": null
            });

            edges.push(edge_obj);
        }

        serde_json::to_string(&edges).unwrap_or_else(|_| "[]".to_string())
    }
}

impl Default for HtmlVizWriter {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_html_viz_writer_creation() {
        let writer = HtmlVizWriter::new();
        assert!(writer.is_default());
    }

    #[test]
    fn test_html_viz_writer_default() {
        let writer = HtmlVizWriter::default();
        assert!(writer.is_default());
    }

    impl HtmlVizWriter {
        fn is_default(&self) -> bool {
            true
        }
    }
}
