#!/bin/bash
sed -i 's/let result = pipeline.run(&input_files)?;/let result = pipeline.run();/g' src/lib.rs
