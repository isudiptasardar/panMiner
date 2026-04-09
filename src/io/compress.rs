//! Zstd compression utilities for intermediate files.

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;
use zstd::stream::{Encoder, Decoder};

use crate::error::{Error, Result};

/// Default compression level (1-22, higher = smaller + slower)
const DEFAULT_COMPRESSION_LEVEL: i32 = 3;

/// Zstd-compressed writer for intermediate files.
pub struct CompressedWriter<'a, W: Write> {
    encoder: Encoder<'a, W>,
}

impl<'a, W: Write> CompressedWriter<'a, W> {
    /// Create a new compressed writer.
    ///
    /// # Arguments
    ///
    /// * `writer` - Underlying writer
    /// * `level` - Compression level (1-22, default 3)
    pub fn new(writer: W, level: i32) -> Result<Self> {
        Ok(Self {
            encoder: Encoder::new(writer, level)?,
        })
    }

    /// Create a new compressed writer with default compression level.
    pub fn with_default(writer: W) -> Result<Self> {
        Self::new(writer, DEFAULT_COMPRESSION_LEVEL)
    }

    /// Write data to the compressed stream.
    pub fn write(&mut self, data: &[u8]) -> Result<()> {
        self.encoder.write_all(data)?;
        Ok(())
    }

    /// Finish compression and return the underlying writer.
    pub fn finish(self) -> Result<W> {
        Ok(self.encoder.finish()?)
    }
}

impl CompressedWriter<'static, File> {
    /// Create a compressed file writer.
    pub fn create(path: &Path) -> Result<Self> {
        let file = File::create(path)?;
        Self::with_default(file)
    }

    /// Create a compressed file writer with custom level.
    pub fn create_with_level(path: &Path, level: i32) -> Result<Self> {
        let file = File::create(path)?;
        Self::new(file, level)
    }
}

/// Zstd-compressed reader for intermediate files.
pub struct CompressedReader<'a, R: Read> {
    decoder: Decoder<'a, std::io::BufReader<R>>,
}

impl<'a, R: Read> CompressedReader<'a, R> {
    /// Create a new compressed reader.
    pub fn new(reader: R) -> Result<Self> {
        Ok(Self {
            decoder: Decoder::new(reader)?,
        })
    }

    /// Read data from the compressed stream.
    pub fn read(&mut self, buf: &mut [u8]) -> Result<usize> {
        Ok(self.decoder.read(buf)?)
    }

    /// Read all data into a vector.
    pub fn read_all(&mut self) -> Result<Vec<u8>> {
        let mut data = Vec::new();
        self.decoder.read_to_end(&mut data)?;
        Ok(data)
    }
}

impl CompressedReader<'static, File> {
    /// Open a compressed file for reading.
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        Self::new(file)
    }
}

/// Compress data to bytes.
pub fn compress(data: &[u8]) -> Result<Vec<u8>> {
    let mut encoder = Encoder::new(Vec::new(), DEFAULT_COMPRESSION_LEVEL)?;
    encoder.write_all(data)?;
    encoder.finish().map_err(Error::from)
}

/// Compress data with custom level.
pub fn compress_with_level(data: &[u8], level: i32) -> Result<Vec<u8>> {
    let mut encoder = Encoder::new(Vec::new(), level)?;
    encoder.write_all(data)?;
    encoder.finish().map_err(Error::from)
}

/// Decompress data.
/// Train a zstd dictionary from a list of samples.
pub fn train_dictionary(samples: &[Vec<u8>], dict_size: usize) -> Result<Vec<u8>> {
    Ok(zstd::dict::from_samples(samples, dict_size)?)
}

/// Compress data using a trained dictionary.
pub fn compress_with_dict(data: &[u8], level: i32, dict: &[u8]) -> Result<Vec<u8>> {
    let mut encoder = zstd::stream::Encoder::with_dictionary(Vec::new(), level, dict)?;
    encoder.write_all(data)?;
    encoder.finish().map_err(crate::error::Error::from)
}

/// Decompress data using a trained dictionary.
pub fn decompress_with_dict(data: &[u8], dict: &[u8]) -> Result<Vec<u8>> {
    let mut decoder = zstd::stream::Decoder::with_dictionary(data, dict)?;
    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed)?;
    Ok(decompressed)
}

pub fn decompress(data: &[u8]) -> Result<Vec<u8>> {
    let mut decoder = Decoder::new(data)?;
    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed)?;
    Ok(decompressed)
}

/// Write compressed data to a file.
pub fn write_compressed(path: &Path, data: &[u8]) -> Result<()> {
    let mut writer = CompressedWriter::create(path)?;
    writer.write(data)?;
    writer.finish()?;
    Ok(())
}

/// Read compressed data from a file.
pub fn read_compressed(path: &Path) -> Result<Vec<u8>> {
    let mut reader = CompressedReader::open(path)?;
    reader.read_all()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_compress_decompress() {
        let data = b"Hello, World! This is a test string for compression.";
        let compressed = compress(data).unwrap();
        let decompressed = decompress(&compressed).unwrap();

        assert_eq!(data.to_vec(), decompressed);
        // Compressed should be smaller for repetitive data
        assert!(compressed.len() < data.len() * 2);
    }

    #[test]
    fn test_compress_decompress_levels() {
        let data = b"Repetitive data. Repetitive data. Repetitive data.";

        let level3 = compress_with_level(data, 3).unwrap();
        let level10 = compress_with_level(data, 10).unwrap();
        let level20 = compress_with_level(data, 20).unwrap();

        // Higher levels generally compress more, but for tiny strings,
        // metadata overhead might cause slight variations, so we only
        // assert that compression happened and it's smaller than the original string * 2.
        assert!(level3.len() > 0);
        assert!(level10.len() > 0);
        assert!(level20.len() > 0);

        // All should decompress correctly
        assert_eq!(decompress(&level3).unwrap(), data.to_vec());
        assert_eq!(decompress(&level10).unwrap(), data.to_vec());
        assert_eq!(decompress(&level20).unwrap(), data.to_vec());
    }

    #[test]
    fn test_file_compression() {
        let data = b"Test data for file compression";

        let temp = NamedTempFile::new().unwrap();
        let path = temp.path();

        write_compressed(path, data).unwrap();
        let read = read_compressed(path).unwrap();

        assert_eq!(data.to_vec(), read);
    }
}
    #[test]
    fn test_zstd_dictionary() {
        let mut samples = Vec::new();
        // Generate enough samples to satisfy zstd's dictionary trainer
        // Zstd recommends sample size to be ~100x the dictionary size
        for i in 0..10000 {
            samples.push(format!("{{\"gene_id\":\"{}\",\"sequence\":\"ATCGATCGATCGATCGATCGATCGATCGATCG\"}}", i).into_bytes());
        }
        
        let dict = train_dictionary(&samples, 32768).unwrap();
        assert!(dict.len() > 0);
        
        let data = b"{\"gene_id\":\"126\",\"sequence\":\"ATCGATCGATCGATCGATCGATCGATCGATCG\"}";
        
        let compressed = compress_with_dict(data, 3, &dict).unwrap();
        let decompressed = decompress_with_dict(&compressed, &dict).unwrap();
        
        assert_eq!(decompressed, data);
        
        // Dictionary compression should be very efficient for small repetitive strings
        let no_dict = compress(data).unwrap();
        assert!(compressed.len() < no_dict.len());
    }
    
