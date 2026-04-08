#!/bin/bash
sed -i 's/zstd = "0.13"/zstd = { version = "0.13", features = ["zdict_builder"] }/g' Cargo.toml
