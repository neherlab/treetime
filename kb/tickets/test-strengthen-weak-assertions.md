# Strengthen weak assertions in collapse, Nexus, and sorted-output tests

## Weak sorted-output assertion

`packages/treetime/src/seq/__tests__/test_mutation.rs:89:`

Uses `windows(2)` loop to check ordering instead of asserting exact expected vector.

## Collapse edge tests assert count not content

`packages/treetime/src/representation/algo/topology_cleanup/__tests__/test_collapse_edge.rs:100,396:`

Success defined as `subs.len()` matching expectation, not exact `Vec<Sub>` content.

## Nexus serialization test uses substring check

`packages/treetime/src/representation/payload/timetree.rs:413:`

Uses `contains()` instead of structured comparison.

## Related issues

- Source: [N-test-quality-deficiencies.md](../issues/N-test-quality-deficiencies.md) -- delete after full resolution
