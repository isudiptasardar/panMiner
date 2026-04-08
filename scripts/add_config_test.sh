#!/bin/bash
cat << 'INNER_EOF' >> src/config.rs

    #[test]
    fn test_compression_level_config() {
        let config = PanminerConfig::new().with_compression_level(10);
        assert_eq!(config.compression_level, 10);

        // Validation should fail for out-of-bounds levels
        let bad_config1 = PanminerConfig::new().with_compression_level(0); // too low (min 1)
        assert!(bad_config1.validate().is_err());

        let bad_config2 = PanminerConfig::new().with_compression_level(23); // too high (max 22)
        assert!(bad_config2.validate().is_err());
    }
INNER_EOF
