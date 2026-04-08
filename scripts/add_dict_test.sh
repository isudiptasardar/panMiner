#!/bin/bash
cat << 'INNER_EOF' >> src/io/compress.rs

    #[test]
    fn test_zstd_dictionary() {
        let samples = vec![
            b"{\"gene_id\":\"123\",\"sequence\":\"ATCGATCG\"}".to_vec(),
            b"{\"gene_id\":\"124\",\"sequence\":\"ATCGATCA\"}".to_vec(),
            b"{\"gene_id\":\"125\",\"sequence\":\"ATCGATCC\"}".to_vec(),
        ];
        
        let dict = train_dictionary(&samples, 1024).unwrap();
        assert!(dict.len() > 0);
        
        let data = b"{\"gene_id\":\"126\",\"sequence\":\"ATCGATCG\"}";
        
        let compressed = compress_with_dict(data, 3, &dict).unwrap();
        let decompressed = decompress_with_dict(&compressed, &dict).unwrap();
        
        assert_eq!(decompressed, data);
        
        // Dictionary compression should be very efficient for small repetitive strings
        let no_dict = compress(data).unwrap();
        assert!(compressed.len() < no_dict.len());
    }
INNER_EOF
