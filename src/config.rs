//! Configuration types for PanMiner.

use std::collections::HashSet;
use std::path::PathBuf;

pub use crate::io::QcMode;

/// Correction mode for error handling.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum CorrectionMode {
    /// Aggressive contamination removal - best for phylogenetic studies
    Strict,
    /// Balanced approach (default)
    #[default]
    Default,
    /// Keep all clusters - useful for rare plasmids
    Sensitive,
}

/// Output format options.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum OutputFormat {
    /// Presence/absence matrix (CSV/TSV)
    Matrix,
    /// Core/accessory alignments (FASTA)
    Alignment,
    /// GML graph format (Cytoscape)
    Graph,
    /// JSON/JSONL format
    Json,
    /// Parquet format (requires --features parquet)
    Parquet,
    /// Interactive HTML visualization (requires --features viz)
    HtmlViz,
    /// Structural variant matrix (CSV)
    Struct,
    /// Structural variant matrix (TSV)
    SVMatrix,
}

impl std::fmt::Display for OutputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OutputFormat::Matrix => write!(f, "matrix"),
            OutputFormat::Alignment => write!(f, "alignment"),
            OutputFormat::Graph => write!(f, "graph"),
            OutputFormat::Json => write!(f, "json"),
            OutputFormat::Parquet => write!(f, "parquet"),
            OutputFormat::HtmlViz => write!(f, "html"),
            OutputFormat::Struct => write!(f, "struct"),
            OutputFormat::SVMatrix => write!(f, "svmatrix"),
        }
    }
}

/// Multiple sequence alignment tool.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum AlignmentTool {
    /// MAFFT alignment (default)
    #[default]
    Mafft,
    /// Clustal Omega
    Clustal,
    /// PRANK
    Prank,
}

/// GPU backend preference.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum GpuBackend {
    /// Auto-detect best available
    #[default]
    Auto,
    /// Force CPU (no GPU)
    Cpu,
    /// Prefer CUDA (NVIDIA)
    Cuda,
    /// Prefer wgpu (cross-platform)
    Wgpu,
}

/// Main configuration for PanMiner pipeline.
#[derive(Clone, Debug)]
pub struct PanminerConfig {
    // Input/output
    /// Input GFF3 files
    pub input_files: Vec<PathBuf>,
    /// Output directory
    pub output_dir: PathBuf,
    /// Temporary directory for intermediate files
    pub temp_dir: PathBuf,

    // Processing
    /// Number of threads (0 = auto-detect)
    pub threads: usize,
    /// Number of genomes to process per chunk (0 = all at once)
    pub chunk_size: usize,
    /// Zstd compression level for intermediate files (1-22, default 3)
    pub compression_level: i32,

    // Clustering
    /// Identity threshold for initial clustering (default: 0.98)
    pub cluster_identity: f32,
    /// Enable MMseqs2 clustering (default: true)
    pub enable_mmseqs: bool,
    /// Path to MMseqs2 binary (None = auto-detect)
    pub mmseqs_path: Option<PathBuf>,
    /// Prefer GPU acceleration when available (default: true)
    pub prefer_gpu: bool,

    // Graph construction
    /// Minimum genome support for a cluster
    pub min_support: usize,

    // Error correction
    /// Correction mode
    pub mode: CorrectionMode,
    /// Contamination removal threshold
    pub contamination_threshold: usize,
    /// Gene family collapse threshold (default: 0.70)
    pub collapse_threshold: f32,

    // Pre-processing QC
    /// Enable pre-processing QC (Mash/CheckM)
    pub enable_qc: bool,
    /// QC mode for filtering
    pub qc_mode: QcMode,
    /// Path to CheckM2 database (optional)
    pub checkm_database_path: Option<PathBuf>,

    // Output
    /// Output formats to generate
    pub outputs: HashSet<OutputFormat>,
    /// Alignment tool to use
    pub alignment_tool: AlignmentTool,
    /// Output file prefix
    pub output_prefix: String,

    // GPU
    /// Force CPU even if GPU available
    pub force_cpu: bool,
    /// Preferred GPU backend
    pub gpu_backend: GpuBackend,

    // Verbosity
    /// Enable verbose logging
    pub verbose: bool,
}

impl Default for PanminerConfig {
    fn default() -> Self {
        Self {
            input_files: Vec::new(),
            output_dir: PathBuf::from("."),
            temp_dir: std::env::temp_dir(),
            threads: 0,
            chunk_size: 100,
            compression_level: 3,
            cluster_identity: 0.98,
            enable_mmseqs: true,
            mmseqs_path: None,
            prefer_gpu: true,
            min_support: 1,
            mode: CorrectionMode::Default,
            contamination_threshold: 2,
            collapse_threshold: 0.70,
            enable_qc: true,
            qc_mode: QcMode::Default,
            checkm_database_path: None,
            outputs: [OutputFormat::Matrix, OutputFormat::Alignment].into_iter().collect(),
            alignment_tool: AlignmentTool::default(),
            output_prefix: String::from("panminer"),
            force_cpu: false,
            gpu_backend: GpuBackend::Auto,
            verbose: false,
        }
    }
}

impl PanminerConfig {
    /// Create a new configuration with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Auto-configure based on system resources.
    pub fn auto_configure() -> Self {
        let threads = std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(num_cpus::get());

        Self {
            threads,
            ..Self::default()
        }
    }

    /// Set input files.
    pub fn with_input_files(mut self, files: Vec<PathBuf>) -> Self {
        self.input_files = files;
        self
    }

    /// Set output directory.
    pub fn with_output_dir(mut self, dir: PathBuf) -> Self {
        self.output_dir = dir;
        self
    }

    /// Set number of threads.
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set chunk size for streaming processing.
    /// Set zstd compression level (1-22).
    pub fn with_compression_level(mut self, level: i32) -> Self {
        self.compression_level = level;
        self
    }

    pub fn with_chunk_size(mut self, size: usize) -> Self {
        self.chunk_size = size;
        self
    }

    /// Set temporary directory for intermediate files.
    pub fn with_temp_dir(mut self, dir: PathBuf) -> Self {
        self.temp_dir = dir;
        self
    }

    /// Set correction mode.
    pub fn with_mode(mut self, mode: CorrectionMode) -> Self {
        self.mode = mode;
        self
    }

    /// Set output formats.
    pub fn with_outputs(mut self, outputs: HashSet<OutputFormat>) -> Self {
        self.outputs = outputs;
        self
    }

    /// Enable or disable MMseqs2 clustering.
    pub fn with_enable_mmseqs(mut self, enable: bool) -> Self {
        self.enable_mmseqs = enable;
        self
    }

    /// Enable or disable GPU acceleration preference.
    pub fn with_prefer_gpu(mut self, prefer: bool) -> Self {
        self.prefer_gpu = prefer;
        self
    }

    /// Force CPU processing.
    pub fn force_cpu(mut self, force: bool) -> Self {
        self.force_cpu = force;
        self
    }

    /// Enable or disable pre-processing QC.
    pub fn with_enable_qc(mut self, enable: bool) -> Self {
        self.enable_qc = enable;
        self
    }

    /// Set QC mode.
    pub fn with_qc_mode(mut self, mode: QcMode) -> Self {
        self.qc_mode = mode;
        self
    }

    /// Set path to CheckM2 database.
    pub fn with_checkm_database_path(mut self, path: PathBuf) -> Self {
        self.checkm_database_path = Some(path);
        self
    }

    /// Validate configuration.
    pub fn validate(&self) -> crate::Result<()> {
        if self.input_files.is_empty() {
            return Err(crate::Error::NoGenomes);
        }

        if self.cluster_identity < 0.5 || self.cluster_identity > 1.0 {
            return Err(crate::Error::Config(format!(
                "cluster_identity must be between 0.5 and 1.0, got {}",
                self.cluster_identity
            )));
        }

        if self.compression_level < 1 || self.compression_level > 22 {
            return Err(crate::Error::Config(format!(
                "compression_level must be between 1 and 22, got {}",
                self.compression_level
            )));
        }

        if self.collapse_threshold < 0.0 || self.collapse_threshold > 1.0 {
            return Err(crate::Error::Config(format!(
                "collapse_threshold must be between 0.0 and 1.0, got {}",
                self.collapse_threshold
            )));
        }

        Ok(())
    }

    /// Get effective thread count.
    pub fn effective_threads(&self) -> usize {
        if self.threads == 0 {
            std::thread::available_parallelism()
                .map(|p| p.get())
                .unwrap_or(1)
        } else {
            self.threads
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = PanminerConfig::default();
        assert_eq!(config.cluster_identity, 0.98);
        assert_eq!(config.chunk_size, 100);
        assert!(config.enable_mmseqs);
    }

    #[test]
    fn test_config_validation() {
        let config = PanminerConfig::default();
        assert!(config.validate().is_err()); // No input files

        let config = config.with_input_files(vec![PathBuf::from("test.gff")]);
        assert!(config.validate().is_ok());
    }

    #[test]
    fn test_config_builder() {
        let config = PanminerConfig::new()
            .with_threads(8)
            .with_chunk_size(50)
            .with_mode(CorrectionMode::Strict);

        assert_eq!(config.threads, 8);
        assert_eq!(config.chunk_size, 50);
        assert_eq!(config.mode, CorrectionMode::Strict);
    }
}
    #[test]
    fn test_compression_level_config() {
        let config = PanminerConfig::new().with_compression_level(10);
        assert_eq!(config.compression_level, 10);

        // Validation should fail for out-of-bounds levels
        let bad_config1 = PanminerConfig::new().with_compression_level(0); // too low (min 1)
        assert!(bad_config1.validate().is_err());

        let bad_config2 = PanminerConfig::new().with_compression_level(23); // too high (max 22)
        assert!(bad_config2.validate().is_err());
    }
