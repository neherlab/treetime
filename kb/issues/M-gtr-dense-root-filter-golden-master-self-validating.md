# Dense GTR golden master bakes v1 root-state filtering into the "v0" oracle

The dense GTR golden master capture script does not capture pure v0 output. It reimplements v0's `nij`/`Ti` accumulation but then applies v1's uninformative-root filter to `root_state`, at [packages/treetime/src/gtr/infer_gtr/**tests**/**fixtures**/gm_infer_gtr_dense_capture](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture) (the `uniform_threshold = 1.0 / n + 1e-10` block). It calls `GTR.infer()` directly rather than `tt.infer_gtr()`, specifically to avoid v0's unfiltered `root_state`.

Consequence: `test_gm_infer_gtr_dense` validates v1 against a v1-derived expectation for the root-state contribution, not against v0. The decision [kb/decisions/gtr-uninformative-root-state-filtering.md](../decisions/gtr-uninformative-root-state-filtering.md) states this openly ("after updating the capture script, all seven real datasets pass at 1e-6"), which means the golden master cannot detect the very divergence the decision introduces. A defect in v1's root filtering would not be caught.

This is an oracle-integrity problem, independent of whether the filtering itself is desirable. The capture script must reflect the reference, not the implementation under test.

## Affected code

- Capture script: [packages/treetime/src/gtr/infer_gtr/**tests**/**fixtures**/gm_infer_gtr_dense_capture](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture)
- Tests: [packages/treetime/src/gtr/infer_gtr/**tests**/test_gm_infer_gtr_dense.rs](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs)
- Decision that introduced it: [kb/decisions/gtr-uninformative-root-state-filtering.md](../decisions/gtr-uninformative-root-state-filtering.md)

## Fix direction

Recapture from true v0 (`tt.infer_gtr()` with unfiltered `root_state`) and compare against v1 run with `filter_uninformative_root = false` for parity, keeping the filtered behavior under an explicit, separately validated path. Requires team decision on the nucleotide filtering policy (the mugration path already defaults to v0 and gates filtering behind `--filter-uninformative-root`).

## Related tickets

- [kb/tickets/test-recapture-dense-gtr-golden-master-from-true-v0.md](../tickets/test-recapture-dense-gtr-golden-master-from-true-v0.md)
