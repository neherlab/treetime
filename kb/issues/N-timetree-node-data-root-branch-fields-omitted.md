# Root node omits placeholder branch-length fields in node data JSON

Augur's refine output writes branch-length fields (`branch_length`, `clock_length`, `mutation_length`) for every node, including the root, where they are placeholder values because the root has no parent edge.

Timetree node data JSON omitted these fields on the root node. The root carried dates and confidence but no branch metrics.

The original reading of the impact was that this had no observed effect, because `augur export v2` reads branch metrics from child edges when accumulating divergence and time and never consumes the root's values. That reasoning was incomplete: augur validates that `mutation_length` is present for every node before anything consumes it, so a missing root field rejects the file outright.

## Status

Resolved by `40796fd8`. `fn build_augur_node_data_json()` now returns `(0.0, Some(0.0), Some(0.0))` for a node with no parent edge instead of `(0.0, None, None)`.

Follow-up: the tests and the doc comment describing the old behaviour were not updated and now fail. See [kb/issues/N-timetree-augur-root-branch-field-tests-stale.md](N-timetree-augur-root-branch-field-tests-stale.md).
