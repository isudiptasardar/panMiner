#!/bin/bash
sed -i '67,72c\
        let mut edges = Vec::new();\
        for e in graph.edges.iter() {\
            let (from, to) = e.key();\
            let genomes = e.value().genomes.clone();\
            edges.push((from.clone(), to.clone(), genomes));\
        }' src/correction/missing.rs
