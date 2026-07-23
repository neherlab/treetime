# Test banner comments replace module structure

Large test files use decorative `// --- ... ---` and `// === ... ===` separators to group unrelated concerns [`packages/treetime/src/ancestral/__tests__/test_fitch_sub.rs#L39`](../../packages/treetime/src/ancestral/__tests__/test_fitch_sub.rs#L39) [`packages/treetime/src/optimize/topology/__tests__/test_prop_merge_shared_mutations.rs#L203`](../../packages/treetime/src/optimize/topology/__tests__/test_prop_merge_shared_mutations.rs#L203).

## Impact

Banner-delimited groups share one import scope and helper namespace even when they test different operations. Test filtering cannot target the visual group, and additions tend to extend already broad files.

## Required structure

Replace each banner only when it corresponds to a real test concern: use a focused file for independently owned behavior or a nested test module when helpers and fixtures are shared. Remove decorative separators that add no ownership information.

## Validation

- Test names and module paths remain descriptive and filterable.
- Shared helpers stay after tests in the project-prescribed `helpers` module.
- Repository search finds no decorative separator comments in maintained Rust tests.
