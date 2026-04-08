#!/bin/bash
sed -i 's/let mmap = MmapFile::open("example.gff")?;/let mmap = MmapFile::open(std::path::Path::new("example.gff")).unwrap();/g' src/io/mmap.rs
