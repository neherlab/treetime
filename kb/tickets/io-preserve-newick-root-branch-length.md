# Preserve root branch lengths in Newick round trips

The Newick parser accepts a root-level branch length but drops it because `struct NewickGraph` has no place to store the value. Preserve this parsed field so reading and writing `(A:0.1,B:0.2):0.5;` does not silently produce `(A:0.1,B:0.2);`.

## Implementation

1. Add `root_branch_length: Option<f64>` to `struct NewickGraph` in [`packages/util-newick/src/types.rs#L7-L16`](../../packages/util-newick/src/types.rs#L7-L16).
2. Store the parsed root branch length in `fn visit_root_branch()` in [`packages/util-newick/src/parse.rs#L78-L107`](../../packages/util-newick/src/parse.rs#L78-L107).
3. Emit the stored value after the root subtree in the Newick writer.
4. Propagate the field through Nexus tree records and every conversion that constructs or consumes `NewickGraph`.
5. Add parser, writer, Newick round-trip, and Nexus round-trip tests for present, absent, zero, and scientific-notation root lengths.

## Acceptance criteria

- Root branch lengths round-trip without value loss.
- Trees without a root branch length retain their existing serialization.
- The new field is preserved by all `NewickGraph` transformations.
- Formatter, linter, and tests pass.

## Related issues

- Source: [kb/issues/N-io-newick-root-branch-length-discarded.md](../issues/N-io-newick-root-branch-length-discarded.md)
