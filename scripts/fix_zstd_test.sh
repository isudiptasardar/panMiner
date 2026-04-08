#!/bin/bash
sed -i 's/for i in 0..100 {/for i in 0..1000 {/g' src/io/compress.rs
