# PanMiner Testing Strategy

## Unit Tests

Unit tests are integrated directly in the source files using Rust's built-in testing framework. Run unit tests with:

```bash
cargo test
```

## Integration Tests

Integration tests are located in the `tests/` directory and test the complete pipeline functionality with sample data.

Run integration tests with:
```bash
cargo test --test integration_test
```

## Development Utilities

Development scripts are located in the `scripts/` directory. These scripts help with:

- Adding new tests (`scripts/add_*_test.sh`)
- Fixing specific issues (`scripts/fix_*.sh`)
- Testing specific components (`scripts/test_pipeline.sh`)

## Manual Testing

To manually test the pipeline with sample data:

```bash
cargo run -- sample1.gff sample2.gff -o output_dir
```

## Continuous Integration

Tests are automatically run on every push to the repository via GitHub Actions.