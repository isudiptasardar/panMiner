#!/bin/bash
sed -i 's/zstd::dict::from_buffers/zstd::dict::from_samples/g' src/io/compress.rs
