# Recapture the dense GTR golden master from true v0

The dense GTR golden master capture script applies v1's uninformative-root filter to `root_state` instead of capturing v0's unfiltered behavior, so `test_gm_infer_gtr_dense` validates v1 against v1. Recapture from the genuine v0 reference and validate the v1 filtering separately.

## Task

- Rewrite [packages/treetime/src/gtr/infer_gtr/**tests**/**fixtures**/gm_infer_gtr_dense_capture](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture) to capture v0's `root_state` unfiltered (via `tt.infer_gtr()` or the exact `treeanc.py:1608-1613` accumulation), with no v1-specific filtering baked in.
- Re-capture `gm_infer_gtr_dense_outputs` and run v1 with `filter_uninformative_root = false` (v0 parity) against it.
- Validate the filtered path (`filter_uninformative_root = true`) separately, with expectations derived from its own specification, not from the v0 oracle.
- Reconcile with the team decision on the nucleotide filtering policy (mugration already defaults to v0 and gates filtering behind `--filter-uninformative-root`).

## Related issues

- Source: [kb/issues/M-gtr-dense-root-filter-golden-master-self-validating.md](../issues/M-gtr-dense-root-filter-golden-master-self-validating.md)
