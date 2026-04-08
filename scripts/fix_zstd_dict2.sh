#!/bin/bash
sed -i 's/train_dictionary(&samples, 1024)/train_dictionary(\&samples, 8192)/g' src/io/compress.rs
