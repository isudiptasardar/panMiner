#!/bin/bash
cat << 'INNER_EOF' > patch.rs
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
INNER_EOF

sed -i '/fn test_zstd_dictionary()/,/^    }/c\
    // Replaced by script\
    ' src/io/compress.rs
sed -i '/\/\/ Replaced by script/r patch.rs' src/io/compress.rs
sed -i '/\/\/ Replaced by script/d' src/io/compress.rs
