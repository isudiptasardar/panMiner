#!/bin/bash
cat << 'INNER_EOF' >> src/pipeline.rs

#[cfg(test)]
mod streaming_tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_dummy_gff(id: &str) -> Result<NamedTempFile> {
        let mut temp = NamedTempFile::with_suffix(".gff").unwrap();
        writeln!(temp, "##gff-version 3").unwrap();
        writeln!(temp, "{}\tProkka\tgene\t100\t200\t.\t+\t.\tID=gene1;product=test", id).unwrap();
        writeln!(temp, "##FASTA").unwrap();
        writeln!(temp, ">{}", id).unwrap();
        writeln!(temp, "ATCGATCGATCGATCG").unwrap();
        Ok(temp)
    }

    #[test]
    fn test_pipeline_chunked_streaming() {
        let gff1 = create_dummy_gff("seq1").unwrap();
        let gff2 = create_dummy_gff("seq2").unwrap();
        let gff3 = create_dummy_gff("seq3").unwrap();

        let temp_dir = tempfile::tempdir().unwrap();
        
        let config = PanminerConfig::default()
            .with_input_files(vec![
                gff1.path().to_path_buf(),
                gff2.path().to_path_buf(),
                gff3.path().to_path_buf(),
            ])
            .with_output_dir(temp_dir.path().to_path_buf())
            .with_chunk_size(2); // 3 files, chunk size 2 -> 2 chunks

        let pipeline = PanminerPipeline::new(config);
        let paths = pipeline.run().expect("Pipeline should run successfully with chunks");

        // The matrix should have 3 genomes
        let matrix_content = std::fs::read_to_string(&paths.matrix).unwrap();
        assert!(matrix_content.contains("seq1"));
        assert!(matrix_content.contains("seq2"));
        assert!(matrix_content.contains("seq3"));
    }
}
INNER_EOF
