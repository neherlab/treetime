# Core tests read dataset files at run time

Some tests of the core library read trees, alignments or metadata from `data/` at run time. The project testing rules forbid file input in unit tests. Most core tests embed their datasets in the test binary through `enum TestDataset` [packages/treetime/src/test_utils/alignment.rs](../../packages/treetime/src/test_utils/alignment.rs), which includes the compressed files with `include_bytes!` and decompresses them in memory. Golden-master tests (`test_gm_*`) read their fixtures by design and are out of scope.

## Current state

These tests still read files:

- `fn test_root_sequence_matches_python_h3n2_na_20()` reads `data/flu/h3n2/20/` with `fn nwk_read_file()` and `fn fasta_read_file()` [packages/treetime/src/ancestral/__tests__/test_python_parity.rs#L51-L62](../../packages/treetime/src/ancestral/__tests__/test_python_parity.rs#L51-L62)
- The four tests of `packages/treetime/src/clock/__tests__/test_clock_dengue100.rs` read the tree and metadata of `data/dengue/100/` through `const DATA_DIR` [packages/treetime/src/clock/__tests__/test_clock_dengue100.rs#L122-L126](../../packages/treetime/src/clock/__tests__/test_clock_dengue100.rs#L122-L126)
- `fn test_contract_dense_sparse_real_gtr()` reads six datasets by path: `flu/h3n2/20`, `ebola/20`, `rsv/a/20`, `dengue/20`, `tb/20` and `lassa/L/50` [packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs#L27-L36](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs#L27-L36)
- `fn test_fitch_gtr_deterministic()` reads five of these datasets by path [packages/treetime/src/gtr/infer_gtr/__tests__/test_fitch_deterministic.rs#L30-L36](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_fitch_deterministic.rs#L30-L36)
- `fn test_nwk_parse_all_datasets()` lists every `tree.nwk` under `data/` and parses each one [packages/treetime/src/io/__tests__/test_nwk.rs#L300-L306](../../packages/treetime/src/io/__tests__/test_nwk.rs#L300-L306)

`enum TestDataset` already embeds `flu/h3n2/20`, `ebola/20` and `rsv/a/20`. The other files are small: the `dengue/100` tests read only `tree.nwk` (3.7 KB) and `metadata.tsv` (55 KB), and the compressed alignments of `dengue/20`, `lassa/L/50` and `tb/20` take 40 KB, 25 KB and 0.6 KB.

## Fix direction

- Add the missing datasets to `enum TestDataset` and read every listed test through it
- `fn test_nwk_parse_all_datasets()` checks that every committed tree parses, so it needs the list of trees. Embed the trees through a list of `include_str!` entries, or move the check to a test that may read files, for example the fixture check of `just fixtures-check`
