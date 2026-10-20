# A test helper builds lineage expectations with non-positive counts

`fn constraint_expectation()` [packages/treetime/src/test_utils/timetree.rs#L26-L43](../../packages/treetime/src/test_utils/timetree.rs#L26-L43) builds the `LineageExpectation` of a tree from its date constraints: only the nodes with a given time get a posterior time, and internal nodes without a date keep the default `NodePosterior`. The expected lineage counts of the result reach zero or fall below it: a `CoalescentPrior` built on such an expectation stops with an internal error at $T_c = 0.1$ and $T_c = 1.0$, while the tests that use the helper (`packages/treetime/src/timetree/__tests__/test_coalescent_timescale.rs`, `packages/treetime/src/coalescent/__tests__/helpers.rs`) pass at the $T_c$ values they choose.

> [!IMPORTANT]
> **Investigation required.** How the undated internal nodes lead to non-positive counts in `LineageExpectation::from_posterior()`, the trees and $T_c$ values of the failing construction, and whether a production path can build an expectation from posteriors that leave internal nodes undated.

## Fix direction

Build the helper's expectation from a complete set of node times (for example the Gaussian timetree of the constraints), so it matches what the pipeline passes to `LineageExpectation`.

## Validation

`CoalescentPrior::new()` on the helper's expectation succeeds for $T_c \in \{0.1, 1, 10\}$ on the trees of the tests above.
