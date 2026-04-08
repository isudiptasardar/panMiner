#!/bin/bash
sed -i 's/\.flat_map(|partial| {/\.flat_map(|partial| {/g' src/graph/concurrent.rs
sed -i 's/partial\.adjacencies\.iter()/partial.adjacencies.par_iter()/g' src/graph/concurrent.rs
sed -i 's/partial\.clusters\.iter()/partial.clusters.par_iter()/g' src/graph/concurrent.rs
