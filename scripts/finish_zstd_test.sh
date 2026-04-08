#!/bin/bash
sed -i 's/let dict = train_dictionary(&samples, 1024).unwrap();/let dict = train_dictionary(\&samples, 1024).unwrap();/g' src/io/compress.rs
