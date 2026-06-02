# Test pipeline error paths

Add unit tests for validation guards in pipeline functions that are currently uncovered.

## Ancestral pipeline

- `site_specific_gtr: true` returns error
- `sample_from_profile != Argmax` with `method != Marginal` returns error
- `method = Joint` returns error listing available methods

Location: `packages/treetime/src/ancestral/pipeline.rs`

## Optimize pipeline

- `damping >= 1.0` returns error
- `damping < 0.0` returns error

Location: `packages/treetime/src/optimize/pipeline.rs`

## Timetree pipeline

- `n_branches_posterior: Some(n)` returns unimplemented error

Location: `packages/treetime/src/timetree/pipeline.rs`

## Related issues

Source: `kb/issues/N-test-coverage-gaps.md`
