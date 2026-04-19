//! MMseqs2-GPU integration for gene clustering.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::{Error, Result};
use crate::graph::{Gene, GeneCluster};
use super::traits::Clusterer;

/// MMseqs2 runner for gene clustering.
///
/// This wraps the MMseqs2 tool for fast sequence clustering,
/// with GPU support when available.
pub struct MMseqsRunner {
    /// Path to MMseqs2 binary
    path: PathBuf,
    /// Whether GPU is available
    use_gpu: bool,
    /// Temporary directory for intermediate files
    tmp_dir: PathBuf,
}

impl MMseqsRunner {
    /// Try to detect MMseqs2 installation.
    ///
    /// Returns `Some(MMseqsRunner)` if MMseqs2 is found in PATH.
    pub fn detect() -> Option<Self> {
        which::which("mmseqs").ok().map(|path| {
            let use_gpu = Self::check_gpu_support(&path);
            Self {
                path,
                use_gpu,
                tmp_dir: std::env::temp_dir(),
            }
        })
    }

    /// Try to detect MMseqs2 at a specific path.
    pub fn from_path(path: PathBuf) -> Result<Self> {
        if !path.exists() {
            return Err(Error::MmseqsNotFound);
        }

        let use_gpu = Self::check_gpu_support(&path);
        Ok(Self {
            path,
            use_gpu,
            tmp_dir: std::env::temp_dir(),
        })
    }

    /// Create a new MMseqs2 runner with explicit settings.
    pub fn new(path: PathBuf, use_gpu: bool, tmp_dir: PathBuf) -> Self {
        Self { path, use_gpu, tmp_dir }
    }

    /// Check if MMseqs2 was compiled with GPU support.
    fn check_gpu_support(path: &Path) -> bool {
        Command::new(path)
            .arg("version")
            .output()
            .map(|out| {
                let stdout = String::from_utf8_lossy(&out.stdout);
                let stderr = String::from_utf8_lossy(&out.stderr);
                stdout.contains("CUDA") || stdout.contains("GPU")
                    || stderr.contains("CUDA") || stderr.contains("GPU")
            })
            .unwrap_or(false)
    }

    /// Check if GPU is available on the system.
    /// Returns true if nvidia-smi is available and reports a GPU.
    pub fn check_system_gpu() -> bool {
        // First check if nvidia-smi is available
        let nvidia_smi = match which::which("nvidia-smi").ok() {
            Some(path) => path,
            None => return false,
        };
        let output = match Command::new(nvidia_smi)
            .arg("--query-gpu=name")
            .arg("--format=csv")
            .output()
        {
            Ok(out) => out,
            Err(_) => return false,
        };

        // Check if any GPU is reported (non-empty output beyond header)
        let stdout = String::from_utf8_lossy(&output.stdout);
        // The CSV output has a header "name" and then GPU names, one per line
        // If there's at least one actual GPU entry (line after header), GPU is available
        let lines: Vec<&str> = stdout.trim().split('\n').collect();
        !lines.is_empty() && (lines.len() > 1 || !lines[0].is_empty())
    }

    /// Check if GPU is available for this runner.
    pub fn has_gpu(&self) -> bool {
        self.use_gpu
    }

    /// Set the temporary directory for intermediate files.
    pub fn set_tmp_dir(&mut self, path: PathBuf) {
        self.tmp_dir = path;
    }

    /// Run MMseqs2 easy-cluster and return the output directory path.
    ///
    /// This is a convenience method that runs the full clustering pipeline.
    /// The actual parsing is done separately to access gene sequences for centroids.
    pub fn easy_cluster(
        &self,
        input: &Path,
        output: &Path,
        identity: f32,
    ) -> Result<PathBuf> {
        let mut cmd = Command::new(&self.path);

        cmd.arg("easy-cluster")
            .arg(input)
            .arg(output)
            .arg(&self.tmp_dir)
            .arg("--min-seq-id")
            .arg(format!("{}", identity))
            .arg("--cluster-mode")
            .arg("2"); // Set cover (greedy)

        if self.use_gpu {
            cmd.arg("--gpu").arg("1");
        }

        let output_result = cmd.output()
            .map_err(|e| Error::Mmseqs(format!("Failed to execute MMseqs2: {}", e)))?;

        if !output_result.status.success() {
            let stderr = String::from_utf8_lossy(&output_result.stderr);
            return Err(Error::Mmseqs(format!("MMseqs2 failed: {}", stderr)));
        }

        Ok(output.to_path_buf())
    }

    /// Parse MMseqs2 cluster output with gene sequences to set centroids.
    fn parse_cluster_output(&self, output_dir: &Path, genes: &[Gene]) -> Result<Vec<GeneCluster>> {
        // MMseqs2 outputs a cluster file with format:
        // representative    member1
        // representative    member2
        // ...
        let cluster_file = output_dir.join("cluster.tsv");

        if !cluster_file.exists() {
            return Err(Error::Mmseqs(format!(
                "Cluster file not found: {:?}",
                cluster_file
            )));
        }

        use std::io::{BufRead, BufReader};
        use std::fs::File;

        // Build gene ID to sequence mapping
        let gene_sequences: std::collections::HashMap<String, Vec<u8>> = genes
            .iter()
            .map(|g| (g.id.to_string(), g.sequence.clone()))
            .collect();

        let file = File::open(&cluster_file)?;
        let reader = BufReader::new(file);

        let mut clusters: std::collections::HashMap<String, GeneCluster> = std::collections::HashMap::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split_whitespace().collect();

            if parts.len() >= 2 {
                let representative = parts[0];
                let member = parts[1];

                let cluster = clusters
                    .entry(representative.to_string())
                    .or_insert_with(|| {
                        let mut c = GeneCluster::new(representative);
                        c.support = 1;
                        c
                    });
                cluster.add_gene(crate::graph::GeneId::new(member));
            }
        }

        // Set centroid sequences from the representative gene
        for cluster in clusters.values_mut() {
            if let Some(seq) = gene_sequences.get(cluster.id.as_str()) {
                cluster.centroids = vec![seq.clone()];
            }
        }

        Ok(clusters.into_values().collect())
    }
}

impl Clusterer for MMseqsRunner {
    fn cluster(&self, genes: &[Gene], identity_threshold: f32) -> Result<Vec<GeneCluster>> {
        // Write genes to FASTA file
        let input_file = self.tmp_dir.join("input.fasta");
        write_genes_to_fasta(genes, &input_file)?;

        // Run clustering
        let output_dir = self.tmp_dir.join("clusters");
        let _clusters = self.easy_cluster(&input_file, &output_dir, identity_threshold)?;

        // Parse cluster output with gene sequences to set centroids
        self.parse_cluster_output(&output_dir, genes)
    }

    fn name(&self) -> &str {
        if self.use_gpu {
            "MMseqs2-GPU"
        } else {
            "MMseqs2-CPU"
        }
    }

    fn is_available(&self) -> bool {
        self.path.exists()
    }
}

/// Write genes to a FASTA file for MMseqs2 input.
fn write_genes_to_fasta(genes: &[Gene], path: &Path) -> Result<()> {
    use std::io::Write;
    use std::fs::File;

    let mut file = File::create(path)?;

    for gene in genes {
        writeln!(file, ">{}", gene.id)?;
        let seq = String::from_utf8_lossy(&gene.sequence);
        writeln!(file, "{}", seq)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_mmseqs() {
        // This test will pass if MMseqs2 is installed
        let runner = MMseqsRunner::detect();
        // Just check it doesn't panic
        if let Some(r) = runner {
            println!("Found MMseqs2 at: {:?}", r.path);
        }
    }
}