# Validate and atomically apply reroot edge splits

Make public reroot topology changes reject invalid fractions and preserve the graph on every error.

## Required changes

- Add a finite `SplitFraction` type bounded to $[0,1]$.
- Use it in split evaluation, root-search results, and graph mutation APIs.
- Resolve every supplied key and required edge/node before the first mutation.
- Replace recoverable `expect` calls with contextual project errors.
- Commit the complete edge split atomically.

## Validation

- Accept $0$, $0.5$, and $1$; reject negative, greater-than-one, NaN, and infinite inputs.
- Assert split lengths are non-negative, finite, and sum to the original length.
- Inject stale node/edge keys and compare the entire graph before and after the error.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-graph-reroot-split-validation-and-atomicity.md](../issues/M-graph-reroot-split-validation-and-atomicity.md)
