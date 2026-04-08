#!/bin/bash
sed -i '/pub chunk_size: usize,/a \    /// Zstd compression level for intermediate files (1-22, default 3)\n    pub compression_level: i32,' src/config.rs
sed -i '/chunk_size: 100,/a \            compression_level: 3,' src/config.rs
sed -i '/pub fn with_chunk_size(mut self, size: usize) -> Self {/i \    /// Set zstd compression level (1-22).\n    pub fn with_compression_level(mut self, level: i32) -> Self {\n        self.compression_level = level;\n        self\n    }\n' src/config.rs
sed -i '/if self.collapse_threshold < 0.0 || self.collapse_threshold > 1.0 {/i \        if self.compression_level < 1 || self.compression_level > 22 {\n            return Err(crate::Error::Config(format!(\n                "compression_level must be between 1 and 22, got {}",\n                self.compression_level\n            )));\n        }\n' src/config.rs
