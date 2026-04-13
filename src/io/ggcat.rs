//! GGCAT colored compacted de Bruijn graph construction wrapper.
//!
//! This module provides data structures and a builder for constructing
//! colored compacted de Bruijn graphs (cDBGs) from genome FASTA files.
//!
//! The `CDBGGraph`, `CDBGStats`, and `GGCATBuilder` types are always available.
//! The actual GGCAT API integration (`build_colored_cdbg`, `compute_stats`)
//! requires the `dbg` feature. Without it, those methods return
//! `Error::FeatureNotEnabled("dbg")`.
//! # Overview
//!
//! - `CDBGGraph`: Result of a cDBG build (paths, k-mer size, unitig/color counts)
//! - `CDBGStats`: Statistics computed from a built cDBG (N50, avg length, color complexity)
//! - `GGCATBuilder`: Builder for configuring and running GGCAT colored cDBG construction
//!
//! # Example
//!
//! ```ignore
//! use panminer::io::GGCATBuilder;
//!
//! let builder = GGCATBuilder::new()
//!     .with_kmer_size(31)
//!     .with_threads(8)
//!     .with_memory_gb(4.0);
//!
//! let genomes = vec![
//!     (PathBuf::from("genome1.fna"), "genome1".to_string()),
//!     (PathBuf::from("genome2.fna"), "genome2".to_string()),
//! ];
//!
//! let graph = builder.build_colored_cdbg(&genomes, 31, Path::new("output/cdbg"))?;
//! let stats = builder.compute_stats(&graph)?;
//! println!("N50: {}, unitigs: {}", stats.n50, stats.num_unitigs);
//! ```

use std::path::{Path, PathBuf};

use crate::error::{Error, Result};

// ---------------------------------------------------------------------------
// CDBGGraph — result of a cDBG build
// ---------------------------------------------------------------------------

/// Result of a GGCAT colored compacted de Bruijn graph build.
///
/// Contains paths to the output graph (GFA format) and color map file,
/// along with metadata about the build parameters and graph size.
#[derive(Debug, Clone)]
pub struct CDBGGraph {
    /// Path to the built graph file (GFA format).
    pub graph_path: PathBuf,
    /// Path to the color map file (mapping unitigs to genome colors).
    pub color_map_path: PathBuf,
    /// Genome identifiers in input order.
    pub color_names: Vec<String>,
    /// K-mer size used for graph construction.
    pub kmer_size: usize,
    /// Number of unitigs in the graph.
    pub num_unitigs: usize,
    /// Number of colors (genomes) in the graph.
    pub num_colors: usize,
}

// ---------------------------------------------------------------------------
// CDBGStats — statistics from a built cDBG
// ---------------------------------------------------------------------------

/// Statistics computed from a built colored compacted de Bruijn graph.
///
/// Includes assembly statistics (N50, unitig counts, total bases) and
/// color-specific metrics (number of colors, average colors per unitig).
#[derive(Debug, Clone, Default)]
pub struct CDBGStats {
    /// N50 of unitig lengths (50% of total bases are in unitigs >= this length).
    pub n50: usize,
    /// Number of unitigs in the graph.
    pub num_unitigs: usize,
    /// Total number of distinct k-mers in the graph.
    pub num_kmers: usize,
    /// Number of colors (genomes) in the graph.
    pub num_colors: usize,
    /// Total bases across all unitigs.
    pub total_bases: usize,
    /// Average unitig length.
    pub avg_unitig_length: f64,
    /// Average number of colors per unitig (color complexity).
    pub color_complexity: f64,
}

// ---------------------------------------------------------------------------
// GGCATBuilder — builder for GGCAT colored cDBG construction
// ---------------------------------------------------------------------------

/// Builder for GGCAT colored compacted de Bruijn graph construction.
///
/// Configures GGCAT runtime parameters (threads, memory, temp directory)
/// and provides methods to build a colored cDBG from genome FASTA files
/// and compute statistics on the resulting graph.
///
/// # Defaults
///
/// - `threads`: 0 (auto-detect)
/// - `memory_gb`: 0.0 (unlimited)
/// - `kmer_size`: 31
/// - `temp_dir`: `None` (system temp)
pub struct GGCATBuilder {
    /// Number of threads for GGCAT (0 = auto-detect).
    threads: usize,
    /// Memory limit in GB for GGCAT (0.0 = unlimited).
    memory_gb: f64,
    /// Default k-mer size for graph construction.
    kmer_size: usize,
    /// Temporary directory for GGCAT intermediate files.
    temp_dir: Option<PathBuf>,
}

impl GGCATBuilder {
    /// Create a new GGCATBuilder with default settings.
    pub fn new() -> Self {
        Self {
            threads: 0,
            memory_gb: 0.0,
            kmer_size: 31,
            temp_dir: None,
        }
    }

    /// Set the number of threads for GGCAT (0 = auto-detect).
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set the memory limit in GB for GGCAT (0.0 = unlimited).
    pub fn with_memory_gb(mut self, memory_gb: f64) -> Self {
        self.memory_gb = memory_gb;
        self
    }

    /// Set the default k-mer size for graph construction.
    pub fn with_kmer_size(mut self, kmer_size: usize) -> Self {
        self.kmer_size = kmer_size;
        self
    }

    /// Set a custom temporary directory for GGCAT intermediate files.
    pub fn with_temp_dir(mut self, path: PathBuf) -> Self {
        self.temp_dir = Some(path);
        self
    }

    /// Build a colored compacted de Bruijn graph from genome FASTA files.
    ///
    /// Takes a list of (FASTA path, genome name) pairs and constructs a
    /// colored cDBG using the GGCAT Rust API. The graph is written to
    /// `output_path` in GFA format with an accompanying color map file.
    ///
    /// # Arguments
    ///
    /// * `genomes` - Slice of (FASTA file path, genome identifier) pairs
    /// * `kmer_size` - K-mer size for graph construction (overrides builder default)
    /// * `output_path` - Output directory for graph and color map files
    ///
    /// # Errors
    ///
    /// Returns `Error::FeatureNotEnabled("dbg")` if the `dbg` feature is not enabled.
    /// Returns `Error::InvalidInput` if the genomes list is empty.
    /// Returns `Error::ggcat_build_failed` if GGCAT fails during construction.
    #[cfg(feature = "dbg")]
    pub fn build_colored_cdbg(
        &self,
        genomes: &[(PathBuf, String)],
        kmer_size: usize,
        output_path: &Path,
    ) -> Result<CDBGGraph> {
        self.build_with_ggcat_api(genomes, kmer_size, output_path)
    }

    /// Build a colored compacted de Bruijn graph.
    ///
    /// Returns `Error::FeatureNotEnabled` when the `dbg` feature is not enabled.
    #[cfg(not(feature = "dbg"))]
    pub fn build_colored_cdbg(
        &self,
        _genomes: &[(PathBuf, String)],
        _kmer_size: usize,
        _output_path: &Path,
    ) -> Result<CDBGGraph> {
        Err(Error::FeatureNotEnabled("dbg".to_string()))
    }

    /// Internal implementation using the GGCAT Rust API.
    ///
    /// Creates a GGCAT instance, builds a colored cDBG from the input
    /// FASTA files, and returns a `CDBGGraph` with paths and metadata.
    #[cfg(feature = "dbg")]
    fn build_with_ggcat_api(
        &self,
        genomes: &[(PathBuf, String)],
        kmer_size: usize,
        output_path: &Path,
    ) -> Result<CDBGGraph> {
        if genomes.is_empty() {
            return Err(Error::InvalidInput(
                "No genomes provided for cDBG construction".to_string(),
            ));
        }

        std::fs::create_dir_all(output_path)?;

        let color_names: Vec<String> = genomes.iter().map(|(_, name)| name.clone()).collect();

        // Build input streams from FASTA files
        let input_streams: Vec<ggcat_api::GeneralSequenceBlockData> = genomes
            .iter()
            .map(|(path, _)| ggcat_api::GeneralSequenceBlockData::FASTA((path.clone(), None)))
            .collect();

        // Configure GGCAT
        let memory = if self.memory_gb > 0.0 {
            self.memory_gb
        } else {
            0.0 // unlimited
        };

        let config = ggcat_api::GGCATConfig {
            temp_dir: self.temp_dir.clone(),
            memory,
            prefer_memory: true,
            total_threads_count: self.threads,
            intermediate_compression_level: None,
            stats_file: None,
            messages_callback: Some(|lvl, msg| match lvl {
                ggcat_api::MessageLevel::Info => tracing::info!("[GGCAT] {}", msg),
                ggcat_api::MessageLevel::Warning => tracing::warn!("[GGCAT] {}", msg),
                ggcat_api::MessageLevel::Error => tracing::error!("[GGCAT] {}", msg),
                ggcat_api::MessageLevel::UnrecoverableError => {
                    tracing::error!("[GGCAT] UNRECOVERABLE: {}", msg)
                }
            }),
        };

        // Create GGCAT instance
        let instance = ggcat_api::GGCATInstance::create(config).map_err(|e| {
            Error::ggcat_build_failed(&format!("Failed to create GGCAT instance: {}", e))
        })?;

        // Build output file path
        let graph_path = output_path.join("cdbg.gfa");

        // Build the colored cDBG
        let built_path = instance
            .build_graph(
                input_streams,
                graph_path.clone(),
                Some(&color_names),
                kmer_size,
                self.threads,
                false,          // forward_only
                None,           // minimizer_length (auto)
                true,           // colors enabled
                1,              // min_multiplicity
                ggcat_api::ExtraElaboration::UnitigLinks,
                None,           // gfa_output_version (default)
            )
            .map_err(|e| Error::ggcat_build_failed(&format!("Graph build failed: {}", e)))?;

        // Get the color map file path
        let color_map_path = ggcat_api::GGCATInstance::get_colormap_file(&built_path);

        // Count unitigs by dumping them
        let num_unitigs = {
            let count = std::sync::atomic::AtomicUsize::new(0);
            instance
                .dump_unitigs(
                    &built_path,
                    kmer_size,
                    None,   // minimizer_length
                    true,   // colors
                    self.threads,
                    false,  // single_thread_output_function
                    |_read, _colors, _same_colors| {
                        count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    },
                )
                .map_err(|e| {
                    Error::ggcat_build_failed(&format!("Unitig dump failed: {}", e))
                })?;
            count.load(std::sync::atomic::Ordering::Relaxed)
        };

        Ok(CDBGGraph {
            graph_path: built_path,
            color_map_path,
            color_names,
            kmer_size,
            num_unitigs,
            num_colors: genomes.len(),
        })
    }

    /// Compute statistics from a built colored cDBG.
    ///
    /// Dumps unitigs from the graph and computes N50, average length,
    /// total bases, and color complexity (average colors per unitig).
    ///
    /// # Errors
    ///
    /// Returns `Error::FeatureNotEnabled("dbg")` if the `dbg` feature is not enabled.
    /// Returns `Error::ggcat_build_failed` if unitig dumping fails.
    #[cfg(feature = "dbg")]
    pub fn compute_stats(&self, graph: &CDBGGraph) -> Result<CDBGStats> {
        use std::sync::atomic::AtomicUsize;

        let lengths = std::sync::Mutex::new(Vec::new());
        let total_color_count = AtomicUsize::new(0);

        let instance = ggcat_api::GGCATInstance::create(ggcat_api::GGCATConfig {
            temp_dir: self.temp_dir.clone(),
            memory: if self.memory_gb > 0.0 {
                self.memory_gb
            } else {
                0.0
            },
            prefer_memory: true,
            total_threads_count: self.threads,
            intermediate_compression_level: None,
            stats_file: None,
            messages_callback: Some(|lvl, msg| match lvl {
                ggcat_api::MessageLevel::Info => tracing::info!("[GGCAT] {}", msg),
                ggcat_api::MessageLevel::Warning => tracing::warn!("[GGCAT] {}", msg),
                ggcat_api::MessageLevel::Error => tracing::error!("[GGCAT] {}", msg),
                ggcat_api::MessageLevel::UnrecoverableError => {
                    tracing::error!("[GGCAT] UNRECOVERABLE: {}", msg)
                }
            }),
        })
        .map_err(|e| {
            Error::ggcat_build_failed(&format!(
                "Failed to create GGCAT instance for stats: {}",
                e
            ))
        })?;

        instance
            .dump_unitigs(
                &graph.graph_path,
                graph.kmer_size,
                None,   // minimizer_length
                true,   // colors
                self.threads,
                false,  // single_thread_output_function
                |read, colors, _same_colors| {
                    lengths.lock().unwrap().push(read.len());
                    total_color_count.fetch_add(
                        colors.len(),
                        std::sync::atomic::Ordering::Relaxed,
                    );
                },
            )
            .map_err(|e| {
                Error::ggcat_build_failed(&format!("Unitig dump for stats failed: {}", e))
            })?;

        let tc = total_color_count.load(std::sync::atomic::Ordering::Relaxed);
        let lv = lengths.into_inner().unwrap();

        compute_cdbg_stats(&lv, tc, graph.num_colors)
    }

    /// Compute statistics from a built colored cDBG.
    ///
    /// Returns `Error::FeatureNotEnabled` when the `dbg` feature is not enabled.
    #[cfg(not(feature = "dbg"))]
    pub fn compute_stats(&self, _graph: &CDBGGraph) -> Result<CDBGStats> {
        Err(Error::FeatureNotEnabled("dbg".to_string()))
    }
}

impl Default for GGCATBuilder {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Stats computation (shared between feature-gated paths and tests)
// ---------------------------------------------------------------------------

/// Compute cDBG statistics from unitig lengths and color counts.
///
/// This is a pure computation function separated from the GGCAT API
/// so it can be tested without the `dbg` feature enabled.
/// Compute cDBG statistics from unitig lengths and color counts.
///
/// This is a pure computation function separated from the GGCAT API
/// so it can be tested without the `dbg` feature enabled. It can also
/// be used directly when unitig data is obtained from other sources.
pub fn compute_cdbg_stats(
    lengths: &[usize],
    total_color_count: usize,
    num_colors: usize,
) -> Result<CDBGStats> {
    if lengths.is_empty() {
        return Ok(CDBGStats {
            num_colors,
            ..Default::default()
        });
    }

    let num_unitigs = lengths.len();
    let total_bases: usize = lengths.iter().sum();
    let avg_unitig_length = total_bases as f64 / num_unitigs as f64;
    let color_complexity = total_color_count as f64 / num_unitigs as f64;

    // Compute N50: sort lengths descending, find length where cumulative sum
    // crosses 50% of total bases
    let mut sorted_lengths: Vec<usize> = lengths.to_vec();
    sorted_lengths.sort_by(|a, b| b.cmp(a)); // descending

    let half_total = total_bases / 2;
    let mut cumulative = 0usize;
    let mut n50 = 0usize;
    for &len in &sorted_lengths {
        cumulative += len;
        if cumulative >= half_total {
            n50 = len;
            break;
        }
    }

    // Estimate num_kmers: each unitig of length L contributes (L - k + 1) k-mers
    // where k is implicit from the graph. We use a simple approximation:
    // total distinct k-mers ≈ sum over unitigs of max(1, L - 30) for k=31
    // This is an approximation; the actual count depends on k-mer overlap.
    let num_kmers = lengths.iter().map(|&l| if l > 31 { l - 31 + 1 } else { 1 }).sum();

    Ok(CDBGStats {
        n50,
        num_unitigs,
        num_kmers,
        num_colors,
        total_bases,
        avg_unitig_length,
        color_complexity,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cdbg_stats_default() {
        let stats = CDBGStats::default();
        assert_eq!(stats.n50, 0);
        assert_eq!(stats.num_unitigs, 0);
        assert_eq!(stats.num_kmers, 0);
        assert_eq!(stats.num_colors, 0);
        assert_eq!(stats.total_bases, 0);
        assert_eq!(stats.avg_unitig_length, 0.0);
        assert_eq!(stats.color_complexity, 0.0);
    }

    #[test]
    fn test_cdbg_graph_fields() {
        let graph = CDBGGraph {
            graph_path: PathBuf::from("/tmp/cdbg.gfa"),
            color_map_path: PathBuf::from("/tmp/cdbg.gfa.colormap"),
            color_names: vec!["genome1".to_string(), "genome2".to_string()],
            kmer_size: 31,
            num_unitigs: 100,
            num_colors: 2,
        };

        assert_eq!(graph.kmer_size, 31);
        assert_eq!(graph.color_names.len(), 2);
        assert_eq!(graph.color_names[0], "genome1");
        assert_eq!(graph.color_names[1], "genome2");
        assert_eq!(graph.num_unitigs, 100);
        assert_eq!(graph.num_colors, 2);
    }

    #[test]
    fn test_ggcat_builder_default() {
        let builder = GGCATBuilder::new();
        assert_eq!(builder.kmer_size, 31);
        assert_eq!(builder.threads, 0);
        assert!((builder.memory_gb - 0.0).abs() < f64::EPSILON);
        assert!(builder.temp_dir.is_none());
    }

    #[test]
    fn test_ggcat_builder_with_methods() {
        let builder = GGCATBuilder::new()
            .with_threads(8)
            .with_memory_gb(4.0)
            .with_kmer_size(51)
            .with_temp_dir(PathBuf::from("/tmp/ggcat"));

        assert_eq!(builder.threads, 8);
        assert!((builder.memory_gb - 4.0).abs() < f64::EPSILON);
        assert_eq!(builder.kmer_size, 51);
        assert_eq!(builder.temp_dir, Some(PathBuf::from("/tmp/ggcat")));
    }

    #[test]
    fn test_build_colored_cdbg_empty_input() {
        let builder = GGCATBuilder::new();
        let result = builder.build_colored_cdbg(&[], 31, Path::new("/tmp/output"));

        // Without dbg feature, it should return FeatureNotEnabled
        #[cfg(not(feature = "dbg"))]
        {
            assert!(result.is_err());
            match result.unwrap_err() {
                Error::FeatureNotEnabled(feat) => assert_eq!(feat, "dbg"),
                _ => panic!("Expected FeatureNotEnabled error"),
            }
        }

        // With dbg feature, it should return InvalidInput
        #[cfg(feature = "dbg")]
        {
            assert!(result.is_err());
            match result.unwrap_err() {
                Error::InvalidInput(msg) => {
                    assert!(msg.contains("No genomes provided"));
                }
                _ => panic!("Expected InvalidInput error"),
            }
        }
    }

    #[test]
    fn test_cdbg_stats_n50_calculation() {
        // Simple case: lengths [100, 50, 30, 20] = total 200
        // Half = 100, cumulative descending: 100 >= 100, so N50 = 100
        let lengths = vec![100, 50, 30, 20];
        let stats = compute_cdbg_stats(&lengths, 8, 4).unwrap();
        assert_eq!(stats.n50, 100);
        assert_eq!(stats.num_unitigs, 4);
        assert_eq!(stats.total_bases, 200);
        assert!((stats.avg_unitig_length - 50.0).abs() < f64::EPSILON);
        assert!((stats.color_complexity - 2.0).abs() < f64::EPSILON); // 8 colors / 4 unitigs

        // Another case: lengths [60, 40, 30, 20, 10] = total 160
        // Half = 80, cumulative descending: 60 < 80, 60+40=100 >= 80, so N50 = 40
        let lengths2 = vec![60, 40, 30, 20, 10];
        let stats2 = compute_cdbg_stats(&lengths2, 10, 5).unwrap();
        assert_eq!(stats2.n50, 40);
        assert_eq!(stats2.total_bases, 160);
        assert!((stats2.avg_unitig_length - 32.0).abs() < f64::EPSILON);
        assert!((stats2.color_complexity - 2.0).abs() < f64::EPSILON); // 10 colors / 5 unitigs
    }

    #[test]
    fn test_compute_cdbg_stats_empty() {
        let stats = compute_cdbg_stats(&[], 0, 3).unwrap();
        assert_eq!(stats.n50, 0);
        assert_eq!(stats.num_unitigs, 0);
        assert_eq!(stats.num_colors, 3);
    }

    #[test]
    fn test_ggcat_builder_default_trait() {
        let builder = GGCATBuilder::default();
        assert_eq!(builder.kmer_size, 31);
        assert_eq!(builder.threads, 0);
    }
}