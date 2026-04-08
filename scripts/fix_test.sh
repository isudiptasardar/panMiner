#!/bin/bash
sed -i 's/let gff1 = create_dummy_gff("seq1").unwrap();/let temp_dir = tempfile::tempdir().unwrap();\n        let dir = temp_dir.path();\n        let gff1 = dir.join("seq1.gff");\n        std::fs::write(\&gff1, "##gff-version 3\nseq1\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene1;product=test\n##FASTA\n>seq1\nATCGATCGATCGATCG\n").unwrap();/' src/pipeline.rs
sed -i 's/let gff2 = create_dummy_gff("seq2").unwrap();/let gff2 = dir.join("seq2.gff");\n        std::fs::write(\&gff2, "##gff-version 3\nseq2\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene1;product=test\n##FASTA\n>seq2\nATCGATCGATCGATCG\n").unwrap();/' src/pipeline.rs
sed -i 's/let gff3 = create_dummy_gff("seq3").unwrap();/let gff3 = dir.join("seq3.gff");\n        std::fs::write(\&gff3, "##gff-version 3\nseq3\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene1;product=test\n##FASTA\n>seq3\nATCGATCGATCGATCG\n").unwrap();/' src/pipeline.rs
sed -i 's/gff1.path().to_path_buf()/gff1.clone()/g' src/pipeline.rs
sed -i 's/gff2.path().to_path_buf()/gff2.clone()/g' src/pipeline.rs
sed -i 's/gff3.path().to_path_buf()/gff3.clone()/g' src/pipeline.rs
sed -i 's/let temp_dir = tempfile::tempdir().unwrap();/\/\/ temp_dir already created above/g' src/pipeline.rs
