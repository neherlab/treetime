# UShER input converts missing branch lengths to explicit zero

## Problem

The MAT contains an embedded Newick tree. Its parser correctly represents branch length as `Option<f64>`, preserving the distinction between an absent length and an explicit zero.

`fn usher_to_graph()` then stores the parsed value in the non-optional `UsherNodeImpl.branch_length: f64` and applies `.unwrap_or(0.0)` in [`packages/treetime-io/src/usher_mat.rs`](../../packages/treetime-io/src/usher_mat.rs). The format adapter therefore converts `None` into `Some(0.0)` when a concrete graph edge is constructed. Exact round trips become impossible, and downstream code can interpret unknown evolutionary distance as observed zero change.

This is an application bug, not a MAT limitation. MAT preserves the embedded Newick spelling, and the common graph edge API already represents branch lengths as `Option<f64>`.

## Required behavior

- Change `UsherNodeImpl.branch_length` to `Option<f64>`.
- Preserve `None`, `Some(0.0)`, and `Some(positive)` without an additional enum. An enum is unnecessary because the application currently has only two semantic states: absent or present; zero is an ordinary present value.
- Preserve every present finite Newick value at the format boundary. Algorithms that interpret branch length as evolutionary distance must apply their existing finite, non-negative precondition before inference; that scientific validation is separate from lossless MAT/Newick parsing.
- Preserve the same optional value when the converter constructs the graph edge.

## Validation

- Round-trip absent, explicit zero, and positive embedded-Newick branch lengths.
- Preserve the same negative-value behavior as the shared Newick parser, and verify that inference boundaries still reject values outside their scientific domain.

## Related tickets

- [kb/tickets/io-preserve-usher-missing-branch-lengths.md](../tickets/io-preserve-usher-missing-branch-lengths.md)
