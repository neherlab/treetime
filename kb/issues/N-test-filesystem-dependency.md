# Unit tests depend on repository-level production datasets

> [!WARNING]
> **Needs review.** The inventory below is incomplete. A search for `data/` paths in Rust test sources finds runtime readers beyond the detailed sections, e.g. `packages/treetime/src/optimize/__tests__/test_convergence_sc2.rs#L28-L29` (`data/sc2/2844`) and the files under "Further candidates". Each candidate still needs classification: unit test to convert, or integration/golden-master test whose contract requires production datasets.

Several unit tests read repository-level `data/` files at runtime when the underlying logic can use test-local fixtures. Their paths use `CARGO_MANIFEST_DIR`, so process working-directory changes do not break them; the defect is coupling unit tests to production dataset availability and layout.

All paths relative to `packages/treetime/src/` unless stated otherwise.

## `ancestral/__tests__/test_python_parity.rs` (1 test)

`test_root_sequence_matches_python_h3n2_na_20` reads `data/flu/h3n2/20/{tree.nwk,aln.fasta.xz}`. The tree is 1.5 KB (one line). The alignment is 20 sequences totaling 28 KB uncompressed FASTA, too large for a const literal but embeddable via `include_str!` with a decompressed fixture file.

## `clock/__tests__/test_clock_dengue100.rs` (2 tests)

Both tests read `data/dengue/100/{tree.nwk,metadata.tsv}`. The tree is 3.7 KB and the metadata is 101 lines (~2 KB). Embeddable as const, but the test's value is tied to this specific biological dataset (force_positive_rate behavior with all 198 root positions negative).

## `optimize/__tests__/test_convergence_sc2.rs` (2 tests)

- `test_convergence_sc2_flu_h3n2_20_converges` reads `data/flu/h3n2/20/{tree.nwk,aln.fasta.xz}`. Same dataset as python parity; shares the conversion approach
- `test_convergence_sc2_sparse_converges_on_sc2_2844` reads `data/sc2/2844/{tree.nwk,aln.fasta.xz}`, a large dataset that cannot become a small test-local fixture

## Further candidates

Rust test sources that reference repository `data/` paths, not yet classified:

- `ancestral/__tests__/test_smoke_gtr_iterations.rs`, `ancestral/__tests__/test_smoke_sample_from_profile.rs`, `optimize/__tests__/test_pipeline_gtr_normalized.rs`, `optimize/__tests__/test_pipeline_reroot.rs`: `data/flu/h3n2/20`
- `gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs`, `gtr/infer_gtr/__tests__/test_fitch_deterministic.rs`: several 20-sample datasets as `rstest` cases
- `ancestral/__tests__/test_derived_mutations.rs`, `timetree/__tests__/test_pipeline.rs` (`data/zika/20`), `io/__tests__/test_nwk.rs` (scans `data/` for `tree.nwk` files)
- Outside `packages/treetime`: `packages/app-commands/src/commands/optimize/__tests__/test_augur_node_data.rs`, `packages/app-commands/src/commands/prune/__tests__/test_mutation_outputs.rs`, `packages/app-napi/src/__tests__/test_backend.rs`, `packages/app-server/src/__tests__/test_routes.rs`, `packages/app-server/src/__tests__/test_confine.rs`, `packages/app-cli/src/__tests__/test_transport_parity.rs`, `packages/app-datasets/src/__tests__/test_examples.rs` (dataset catalog, whose contract is the `data/` tree), and further `packages/app-commands` tests that name `data/` paths, possibly only as strings

## Fix

Move small biological inputs into test-local committed fixtures loaded at compile time, or reclassify tests whose contract requires production datasets as integration or golden-master tests.

## Acceptance criteria

- The inventory covers every unit-test source in the repository, not only the examples above
- No unit test reads repository-level production `data/` at runtime
- Small text inputs use test-local fixtures; unit tests do not need compressed production datasets at runtime
- Shared fixtures are stored once and reused by every test of the same dataset
- The biological scenarios and asserted properties stay unchanged
- Unit tests pass with the production `data/` tree temporarily unavailable
