# Dense GTR golden master compares only the filtered root state, not v0 output

The dense GTR golden master capture script does not capture pure v0 output. It reimplements v0's `nij`/`Ti` accumulation, then applies v1's uninformative-root filter to `root_state` [`packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture#L160-L166`](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture#L160-L166): a root position counts only when the maximum of its profile is above `1/n + 1e-10`. It calls `GTR.infer()` directly [`packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture#L174-L179`](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture#L174-L179) instead of `tt.infer_gtr()`, to avoid v0's unfiltered `root_state`.

The Python filter is an independent reimplementation of the Rust filter (`fn is_profile_informative()` [packages/treetime/src/gtr/infer_gtr/common.rs#L245-L248](../../packages/treetime/src/gtr/infer_gtr/common.rs#L245-L248), applied in [packages/treetime/src/partition/marginal/shared/data.rs#L40-L48](../../packages/treetime/src/partition/marginal/shared/data.rs#L40-L48)), so `test_gm_infer_gtr_dense` catches a defect in v1's filter code. It does not check v0 parity: no test compares v1 dense GTR inference with the unfiltered `root_state` that v0 computes. The decision [kb/decisions/gtr-uninformative-root-state-filtering.md](../decisions/gtr-uninformative-root-state-filtering.md) records the size of that divergence only from a measurement made before the capture script was changed (about 29% shift in W matrix elements on `lassa_L_50`); after the change all real datasets pass at 1e-6, so the golden master no longer shows the divergence the decision introduces.

Dense partitions always filter: `DenseInputs` is built with `filter_uninformative_root: true` [packages/treetime/src/partition/marginal/dense/partition.rs#L56-L59](../../packages/treetime/src/partition/marginal/dense/partition.rs#L56-L59), and no flag changes it. Mugration exposes the filter as `--filter-uninformative-root`, off by default (v0 behavior) [packages/app-commands/src/commands/mugration/args.rs#L336-L342](../../packages/app-commands/src/commands/mugration/args.rs#L336-L342).

## Affected code

- Capture script: [`packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture`](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture)
- Tests: [`packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs)
- Decision that introduced the filter: [kb/decisions/gtr-uninformative-root-state-filtering.md](../decisions/gtr-uninformative-root-state-filtering.md)

## Fix direction

Add an unfiltered v0-parity case: capture from true v0 (`tt.infer_gtr()` with unfiltered `root_state`) and compare with v1 counts computed with `filter_uninformative_root = false`. Keep the current filtered capture as the test of the filtered path.

## Proposed implementation plan

The capture script, `fn is_profile_informative()`, the root-state loop in `shared/data.rs`, the hard-coded `filter_uninformative_root: true` in `dense/partition.rs`, and the mugration flag are at the cited lines.

1. Extend `gm_infer_gtr_dense_capture` with a second output per case: call v0 `tt.infer_gtr(pc=1.0, ...)` directly, so v0 computes its own unfiltered `root_state`, and store the inferred `W`, `pi`, `mu` under a separate key (for example `unfiltered`) in `gm_infer_gtr_dense_outputs.json`. Keep the current filtered output as it is. Import v0 code only; do not reimplement the accumulation for this key.
2. Add a test-only way to build the dense inputs with `filter_uninformative_root: false` (for example a constructor argument or a builder on `DenseInputs` reachable from `__tests__`), without a CLI flag, because dense partitions always filter in production.
3. In `test_gm_infer_gtr_dense.rs`, add `test_gm_infer_gtr_dense_unfiltered_synthetic` and `..._real` (rstest cases equal to the existing ones) that run `fn setup_dense_partition()` with the filter off and compare `W`, `pi`, `mu` with the `unfiltered` outputs at 1e-6 (project limit). Oracle: v0 `TreeAnc.infer_gtr()` output, captured by the script.
4. Rerun the capture in the container (v0 Python CLI per `cli-verification.instructions.md`) and review the new outputs before commit.
5. Update [kb/decisions/gtr-uninformative-root-state-filtering.md](../decisions/gtr-uninformative-root-state-filtering.md) only with the measured divergence of filtered vs unfiltered output, with user approval (decision files change only with consent).

Validation: `just test-rs -E 'package(treetime)'`, then `just test-rs`; see red by setting the unfiltered inputs to filter and confirming the new test fails on a dataset with ambiguous root columns (for example `lassa_L_50`).
