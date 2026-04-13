# Sparse GTR test oracle is internally consistent but not independently checked

## `test_get_mutation_counts_sparse` uses a shared-code oracle

`test_get_mutation_counts_sparse` at [packages/treetime/src/gtr/infer_gtr/**tests**/test_sparse.rs#L98](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_sparse.rs#L98) asserts that `get_mutation_counts_sparse` matches counts derived from `expected_mutation_counts_from_reconstructed_sequences` at [packages/treetime/src/gtr/infer_gtr/**tests**/test_sparse.rs#L31](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_sparse.rs#L31). The oracle calls `ancestral_reconstruction_marginal` on the same partition passed to the mutation-count function under test.

Both paths share the canonical-state reconstruction primitives. A defect in the shared code (mask precedence, ambiguous-leaf resolution, canonical-state walk) would shift both `counts_actual` and `counts_expected` by the same amount and the assertion would still pass.

The test verifies internal consistency. It does not verify correctness against an external oracle.

## `test_infer_gtr_sparse` does not call `infer_gtr_sparse`

`test_infer_gtr_sparse` at [packages/treetime/src/gtr/infer_gtr/**tests**/test_sparse.rs#L143](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_sparse.rs#L143) constructs both `actual` and `expected` by calling `infer_gtr_impl(...)` at [packages/treetime/src/gtr/infer_gtr/**tests**/test_sparse.rs#L175](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_sparse.rs#L175). It never calls the public entry point `infer_gtr_sparse` at [packages/treetime/src/gtr/infer_gtr/sparse.rs#L13](../../packages/treetime/src/gtr/infer_gtr/sparse.rs#L13).

The wrapper's default-option wiring and `GTR::new` assembly are not exercised. The body effectively reduces to the first test's claim because `infer_gtr_impl` is deterministic.

## Oracle candidates for replacement

The sparse tests should compare against an oracle that does not share sparse reconstruction code. Options:

- **O1. Dense inference on the same alignment**. `get_mutation_counts_dense` exists at [packages/treetime/src/gtr/infer_gtr/dense.rs#L136](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L136). Build a dense partition from the same compressed data, run the dense mutation counter, compare. The dense path uses full per-site profiles and does not share the sparse state-walk primitives; an equivalence assertion within floating-point tolerance is meaningful. Existing dense/sparse equivalence pattern at [packages/treetime/src/commands/optimize/**tests**/test_dense_sparse_equivalence/test_dense_sparse_equivalence_initial.rs#L70](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_initial.rs#L70) is the closest template.
- **O2. Hand-computed counts for a small fixture**. Use a fixture with 3-5 sequences of length 8-16 and derive expected `n_ij` and `T_i` by hand from the known topology and sequences. Lower fixture size keeps the hand calculation tractable. Suitable when the sparse-vs-dense mapping has its own equivalence mismatches that would undermine O1.
- **O3. v0 Python reference**. Run the equivalent v0 flow (`./dev/docker/python treetime ancestral` on the same fixture, extract mutations from the Python reconstruction) and capture counts as golden values. This is the canonical port oracle. Requires a small capture script following the v0 bootstrap pattern in project rules (seed RNG, fixed alphabet).

For `test_infer_gtr_sparse` specifically, the wrapper-level assertion should call `infer_gtr_sparse(graph, partition, options)` directly and compare the returned `GTR` against either the result of running dense inference on the same input (O1) or a v0 reference capture (O3). Avoid deriving expected values by re-running `infer_gtr_impl` on sparse-derived counts, which is what the current test does.

## Proposed tests

- Replace the reconstructed-sequence oracle in `test_get_mutation_counts_sparse` with O1 (dense equivalence) or O3 (v0 capture).
- Rewrite `test_infer_gtr_sparse` to call `infer_gtr_sparse` directly and compare `W`, `pi`, `mu` to the chosen oracle.

Fixture sizing recommendation: 3-5 sequences, length 16-32, no ambiguity except the ones explicitly being tested. Keep ambiguous-site fixtures separate from the main assertion to isolate what is being checked.

## What it catches

- Drift in `get_mutation_counts_sparse` hidden by a shared reconstruction defect
- Breakage in the `infer_gtr_sparse` wrapper's option defaults, input wiring, or output assembly
- Divergence between sparse and dense mutation counting contracts (O1 path)
- Regression against v0 reference output (O3 path)
