#!/bin/bash
sed -i '170,183c\
        let mut neighbors_a = Vec::new();\
        for entry in graph.edges.iter() {\
            let (from, to) = entry.key();\
            if from == id_a {\
                neighbors_a.push(to.clone());\
            } else if to == id_a {\
                neighbors_a.push(from.clone());\
            }\
        }' src/correction/fragment.rs
