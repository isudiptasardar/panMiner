#!/bin/bash
sed -i '/fn test_zstd_dictionary()/,/^    }/c\
    #[test]\
    fn test_zstd_dictionary() {\
        let mut samples = Vec::new();\
        for i in 0..1000 {\
            samples.push(format!("{{\"gene_id\":\"{}\",\"sequence\":\"ATCGATCGATCGATCGATCGATCGATCGATCG\"}}", i).into_bytes());\
        }\
        let dict = train_dictionary(\&samples, 8192).unwrap();\
        assert!(dict.len() > 0);\
        \
        let data = b"{\"gene_id\":\"126\",\"sequence\":\"ATCGATCGATCGATCGATCGATCGATCGATCG\"}";\
        let compressed = compress_with_dict(data, 3, \&dict).unwrap();\
        let decompressed = decompress_with_dict(\&compressed, \&dict).unwrap();\
        \
        assert_eq!(decompressed, data);\
    }' src/io/compress.rs
