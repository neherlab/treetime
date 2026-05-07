# rayon::ThreadPoolBuilder::build_global called per-test

Calling `build_global` multiple times causes nondeterministic test failures because the global pool can only be initialized once.

v1: [`packages/treetime/src/graph/__tests__/graph.rs#L148`](../../packages/treetime/src/graph/__tests__/graph.rs#L148), [`graph.rs#L167`](../../packages/treetime/src/graph/__tests__/graph.rs#L167) and multiple ancestral test files.

## Related issues

- Source: `do../issues/N-code-quality-conventions.md` -- delete after full resolution
