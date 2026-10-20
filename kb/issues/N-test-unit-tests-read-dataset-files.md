# Unit tests read dataset files at run time

Some unit tests read trees, alignments, or metadata from the repository `data/` folder at run time. The project testing rules forbid file input in unit tests, and the tests break when the dataset layout changes. Most core tests embed their datasets in the test binary through `enum TestDataset` [packages/treetime/src/test_utils/alignment.rs](../../packages/treetime/src/test_utils/alignment.rs), which includes the compressed files with `include_bytes!` and decompresses them in memory. Golden-master tests (`test_gm_*`) read their fixtures by design and are out of scope.

## Core tests that still read files

- `fn test_root_sequence_matches_python_h3n2_na_20()` reads `data/flu/h3n2/20/` with `fn nwk_read_file()` and `fn fasta_read_file()` [packages/treetime/src/ancestral/__tests__/test_python_parity.rs#L51-L62](../../packages/treetime/src/ancestral/__tests__/test_python_parity.rs#L51-L62)
- The four tests of `packages/treetime/src/clock/__tests__/test_clock_dengue100.rs` read the tree and metadata of `data/dengue/100/` through `const DATA_DIR` [packages/treetime/src/clock/__tests__/test_clock_dengue100.rs#L118](../../packages/treetime/src/clock/__tests__/test_clock_dengue100.rs#L118). The tests check the behavior of `force_positive_rate` on this dataset, where every root position gives a negative rate, so the scenario must stay the same
- `fn test_contract_dense_sparse_real_gtr()` reads six datasets by path: `flu/h3n2/20`, `ebola/20`, `rsv/a/20`, `dengue/20`, `tb/20` and `lassa/L/50` [packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs#L27-L36](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs#L27-L36)
- `fn test_fitch_gtr_deterministic()` reads five of these datasets by path [packages/treetime/src/gtr/infer_gtr/__tests__/test_fitch_deterministic.rs#L30-L36](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_fitch_deterministic.rs#L30-L36)
- `fn test_nwk_parse_all_datasets()` lists every `tree.nwk` under `data/` and parses each one [packages/treetime/src/io/__tests__/test_nwk.rs#L300-L306](../../packages/treetime/src/io/__tests__/test_nwk.rs#L300-L306)

`enum TestDataset` already embeds `flu/h3n2/20`, `ebola/20`, `rsv/a/20`, `zika/20`, `zika/86` and `sc2/2844`. The other files are small: the `dengue/100` tests read only `tree.nwk` (3.7 KB) and `metadata.tsv` (55 KB), and the compressed alignments of `dengue/20`, `lassa/L/50` and `tb/20` take 40 KB, 25 KB and 0.6 KB.

## Tests of the application crates

Tests of `app-commands`, `app-backend`, `app-server`, `app-napi`, and `app-cli` also name `data/` paths. Some of them read the files, for example `packages/app-commands/src/commands/optimize/__tests__/test_augur_node_data.rs`, `packages/app-commands/src/commands/prune/__tests__/test_mutation_outputs.rs`, and `packages/app-commands/src/commands/timetree/__tests__/test_pipeline.rs`; others only use the path as a string. Each one needs a classification: a unit test to convert, or an integration test whose contract requires the example datasets, such as the dataset catalog test `packages/app-datasets/src/__tests__/test_examples.rs`.

## Fix direction

- Add the missing datasets to `enum TestDataset` and read every listed core test through it
- `fn test_nwk_parse_all_datasets()` checks that every committed tree parses, so it needs the list of trees. Embed the trees through a list of `include_str!` entries, or move the check to a test that may read files, for example the fixture check of `just fixtures-check`
- Classify the tests of the application crates and convert the unit tests among them

## Acceptance criteria

- No unit test reads the repository `data/` folder at run time
- Shared fixtures are stored once and reused by every test of the same dataset
- The biological scenarios and asserted properties stay unchanged
- Unit tests pass with the `data/` folder temporarily unavailable
