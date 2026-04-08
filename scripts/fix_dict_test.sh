#!/bin/bash
sed -i 's/train_dictionary(&samples, 1024)/train_dictionary(\&samples, 64)/g' src/io/compress.rs
