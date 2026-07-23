# Mugration analysis interface exposes policy wiring

`fn execute_mugration()` accepts eleven positional arguments spanning traits, missing-data policy, weighting, regularization, stopping, bias correction, and behavior switches [`packages/treetime/src/mugration/mugration.rs#L80`](../../packages/treetime/src/mugration/mugration.rs#L80). Tests repeat the positional ordering, and `mugration::pipeline` only re-exports the function as `run` [`packages/treetime/src/mugration/pipeline.rs#L1`](../../packages/treetime/src/mugration/pipeline.rs#L1).

The public interface mirrors implementation wiring instead of representing a complete analysis request. Define a typed mugration request whose states make policy combinations legible, and keep file parsing and output adapters outside the analysis boundary.

## Parameter groups

- Observations: graph, trait map, and attribute name.
- Missing-data policy: missing marker and minimum weight coverage.
- Model initialization: optional weights and pseudo-count.
- Refinement: iteration count and sampling-bias correction.
- Reconstruction policy: initial-frequency smoothing and uninformative-root filtering.

Several adjacent arguments share primitive types, and tests repeat tails such as optional scalars, integers, and booleans whose meaning is not visible at the call site. The alias-only `pipeline` module adds another name without hiding any of this sequence.

## Required module boundary

Make mugration own a validated analysis request, state discovery, equilibrium initialization, GTR refinement, reconstruction, and result assembly. Keep graph/data inputs distinct from inference policy where that improves reuse. Remove the alias-only pipeline module once callers use the cohesive operation directly.

File parsing, CLI defaults, and output serialization remain adapter responsibilities.

## Validation

- Construction tests reject invalid thresholds, iteration values, and inconsistent weighting configuration.
- Golden-master and unit tests preserve missing-state handling, inferred frequencies, rates, reconstruction, and bias correction.
- Callers no longer pass positional policy scalars or booleans.
- Direct analysis tests use in-memory traits and results; command integration tests cover I/O separately.

## Related issues

- [H-core-command-module-shared-ops-entanglement.md](H-core-command-module-shared-ops-entanglement.md)
