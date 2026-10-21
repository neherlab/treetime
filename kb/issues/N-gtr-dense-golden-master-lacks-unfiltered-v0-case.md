# Dense GTR golden master compares only the filtered root state, not v0 output

The dense GTR golden master capture script does not capture pure v0 output. It reimplements v0's `nij`/`Ti` accumulation, then applies v1's uninformative-root filter to `root_state` [packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture#L158-L166](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture#L158-L166): a root position counts only when the maximum of its profile is above `1/n + 1e-10`. It calls `GTR.infer()` directly [packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture#L174-L179](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture#L174-L179) instead of `tt.infer_gtr()`, to avoid v0's unfiltered `root_state`.

The Python filter is an independent reimplementation of the Rust filter (`fn is_profile_informative()` [packages/treetime/src/gtr/infer_gtr/common.rs#L127-L130](../../packages/treetime/src/gtr/infer_gtr/common.rs#L127-L130), applied in [packages/treetime/src/partition/marginal/shared/data.rs#L45-L51](../../packages/treetime/src/partition/marginal/shared/data.rs#L45-L51)), so `test_gm_infer_gtr_dense` catches a defect in v1's filter code. It does not check v0 parity: no test compares v1 dense GTR inference with the unfiltered `root_state` that v0 computes. The decision [kb/decisions/gtr-uninformative-root-state-filtering.md](../decisions/gtr-uninformative-root-state-filtering.md) records the size of that divergence only from a measurement made before the capture script was changed (about 29% shift in W matrix elements on `lassa_L_50`); after the change all real datasets pass at 1e-6, so the golden master no longer shows the divergence the decision introduces.

Dense partitions always filter: `DenseInputs` is built with `filter_uninformative_root: true` [packages/treetime/src/partition/marginal/dense/partition.rs#L59-L62](../../packages/treetime/src/partition/marginal/dense/partition.rs#L59-L62), and no flag changes it. Mugration exposes the filter as `--filter-uninformative-root`, off by default (v0 behavior) [packages/app-commands/src/commands/mugration/args.rs#L170-L176](../../packages/app-commands/src/commands/mugration/args.rs#L170-L176).

## Affected code

- Capture script: [packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture)
- Tests: [packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs)
- Decision that introduced the filter: [kb/decisions/gtr-uninformative-root-state-filtering.md](../decisions/gtr-uninformative-root-state-filtering.md)

## Fix direction

Add an unfiltered v0-parity case: capture from true v0 (`tt.infer_gtr()` with unfiltered `root_state`) and compare with v1 counts computed with `filter_uninformative_root = false`. Keep the current filtered capture as the test of the filtered path.
