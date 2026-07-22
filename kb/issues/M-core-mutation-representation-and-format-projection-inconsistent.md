# Mutation representation and format projection are inconsistent

Mutation coordinates, string rendering, and tree-format projection do not share one application-wide contract. The same biological mutation can change position base, textual spelling, or representable variants depending on the code path.

> [!NOTE]
> The tree-output refactor removed the `tree_ir` layer and now writes each format directly from the graph in [`packages/treetime/src/commands/shared/tree_output.rs`](../../packages/treetime/src/commands/shared/tree_output.rs). This relocated the format-projection sites but did not, by itself, establish a single shared mutation contract. Which of the sub-points below the refactor already resolved is **not yet confirmed**; re-audit the writers in `tree_output.rs` against each point.

## Problem

- Core and sparse mutation structures mix zero-based and one-based positions.
- Substitution strings are assembled at multiple output sites instead of by the mutation value itself.
- Insertions and deletions reconstructed by ancestral inference may not reach every format's output.
- UShER conversion narrows positions through `usize`/`i32` casts and reconstructs one-based format positions ad hoc.

Affected boundaries include `enum Mutation` and sparse mutation types in [`packages/treetime/src/partition`](../../packages/treetime/src/partition), and the Auspice, PhyloXML, and UShER writers in [`packages/treetime/src/commands/shared/tree_output.rs`](../../packages/treetime/src/commands/shared/tree_output.rs).

## Required invariant

Let $p_0$ be the internal zero-based position and $p_1$ the one-based position required by an external format. Conversion occurs only at the format boundary:

$$p_1 = p_0 + 1$$

The inverse conversion is defined only for $p_1 \ge 1$:

$$p_0 = p_1 - 1$$

Every narrowing conversion must also prove the target integer range. Position zero, negative positions, and overflow return contextual errors.

## Design axes

### Coordinate storage

- O1. Store zero-based positions throughout the application, including TreeIR. This matches sequence/array indexing and makes one-based values an external serialization concern.
- O2. Store format-native positions in each representation. This preserves wire spelling but spreads conversion rules through algorithms.

Recommendation: O1. Internal positions are zero-based throughout the application; one-based positions exist only while parsing or serializing a format that requires them.

### Mutation rendering

- O1. Put canonical mutation-string methods on the mutation type and give substitutions, insertions, and deletions explicit variants.
- O2. Keep per-writer formatting helpers. This permits format-specific control but duplicates shared spelling and coordinate conversion.

Recommendation: O1 for the canonical application string, with checked format-specific serializers calling the same typed mutation API.

### Format projection

- O1. Project substitutions, insertions, and deletions only where the external format has a verified lossless representation.
- O2. Define a documented lossy projection for a specific format. This changes scientific output and requires explicit approval.

Recommendation: O1 for verified mappings. The UShER writer now rejects unrepresentable (amino-acid, indel) mutations with an explicit error rather than dropping them silently. The PhyloXML mutation-property loss policy remains open in [N-io-phyloxml-mutation-property-contract-undecided.md](N-io-phyloxml-mutation-property-contract-undecided.md) until its mapping is approved. Silent omission is invalid.

## Recommendation

Use zero-based typed mutations throughout the application, put canonical rendering on the mutation type, and perform checked conversion only at format boundaries. Implement only lossless, verified format mappings; unresolved UShER and PhyloXML loss policies remain separate decisions.

## Validation

- Unit tests for substitution, insertion, and deletion strings at internal positions $0$, $1$, and a large valid coordinate.
- Property tests for checked zero-based/one-based round trips.
- End-to-end projection tests for every verified mutation/format mapping.
- Explicit rejection tests for zero, negative, overflow, and format-unrepresentable mutations.

## Related issues

- [N-io-phyloxml-mutation-property-contract-undecided.md](N-io-phyloxml-mutation-property-contract-undecided.md)
