#!/bin/bash
sed -i 's/let parser = GffParser::open("genome.gff",/let parser = GffParser::open(std::path::Path::new("genome.gff"),/' src/io/gff.rs
sed -i 's/let parser = FastaParser::open("sequences.fasta")?;/let parser = FastaParser::open(std::path::Path::new("sequences.fasta")).unwrap();/' src/io/fasta.rs
sed -i 's/let mmap = MmapFile::open("large_file.txt")?;/let mmap = MmapFile::open(std::path::Path::new("large_file.txt")).unwrap();/' src/io/mmap.rs
sed -i 's/let bytes = mmap.as_bytes();/let bytes = mmap.as_bytes();/' src/io/mmap.rs
sed -i 's/let genes = parser.parse_genes()?;/let genes = parser.parse_genes().unwrap();/' src/io/gff.rs
sed -i 's/let parser = GffParser::open(std::path::Path::new("genome.gff"), GenomeId::new("sample1"))?;/let parser = GffParser::open(std::path::Path::new("genome.gff"), GenomeId::new("sample1")).unwrap();/' src/io/gff.rs
sed -i 's/let sequences = parser.parse_all()?;/let sequences = parser.parse_all().unwrap();/' src/io/fasta.rs
