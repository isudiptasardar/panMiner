//! Memory-mapped file handling for zero-copy parsing.

use std::fs::File;
use std::path::Path;
use memmap2::Mmap;

use crate::error::Result;

/// A memory-mapped file for efficient zero-copy parsing.
///
/// This wraps `memmap2::Mmap` and provides safe access to file contents
/// without loading the entire file into memory.
pub struct MmapFile {
    mmap: Mmap,
    _file: File,
}

impl MmapFile {
    /// Open a file and memory-map it.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the file to memory-map
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or mapped.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use panminer::io::MmapFile;
    ///
    /// # fn main() -> panminer::error::Result<()> {
    /// let mmap = MmapFile::open(std::path::Path::new("example.gff"))?;
    /// let bytes = mmap.as_bytes();
    /// # Ok(())
    /// # }
    /// ```
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(Self { mmap, _file: file })
    }

    /// Get the entire file contents as a byte slice.
    pub fn as_bytes(&self) -> &[u8] {
        &self.mmap
    }

    /// Get the length of the file in bytes.
    pub fn len(&self) -> usize {
        self.mmap.len()
    }

    /// Check if the file is empty.
    pub fn is_empty(&self) -> bool {
        self.mmap.is_empty()
    }

    /// Iterate over lines in the file without allocating strings.
    ///
    /// This returns an iterator over byte slices, each representing a line
    /// without the newline character. Handles both Unix (\n) and Windows (\r\n) endings.
    pub fn lines(&self) -> impl Iterator<Item = &[u8]> {
        self.mmap.split(|&b| b == b'\n')
            .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
    }

    /// Get a line by line number (0-indexed).
    ///
    /// Returns `None` if the line number is out of bounds.
    pub fn get_line(&self, line_num: usize) -> Option<&[u8]> {
        self.lines().nth(line_num)
    }

    /// Count the number of lines in the file.
    pub fn count_lines(&self) -> usize {
        self.lines().count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_mmap_file() -> Result<()> {
        let mut temp = NamedTempFile::new()?;
        write!(temp, "line 1\nline 2\nline 3")?;

        let mmap = MmapFile::open(temp.path())?;
        assert_eq!(mmap.len(), 20);
        assert!(!mmap.is_empty());

        Ok(())
    }

    #[test]
    fn test_lines() -> Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "line 1\nline 2\nline 3")?;

        let mmap = MmapFile::open(temp.path())?;
        let lines: Vec<&[u8]> = mmap.lines().filter(|l| !l.is_empty()).collect();
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], b"line 1");
        assert_eq!(lines[2], b"line 3");

        Ok(())
    }
}