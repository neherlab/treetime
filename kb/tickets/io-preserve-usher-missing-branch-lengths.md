# Preserve missing UShER branch lengths

Retain the distinction between an absent branch length and an explicit zero in the embedded MAT Newick tree.

## Required changes

- Change `UsherNodeImpl.branch_length` from `f64` to `Option<f64>`.
- Remove the `.unwrap_or(0.0)` conversion in `fn usher_to_graph()`.
- Pass the optional value unchanged through every `UsherRead` implementation into `HasBranchLength::set_branch_length()` or the edge constructor.
- Preserve the shared Newick parser's present-value behavior; do not introduce MAT-only branch-length validation.

## Validation

- Round-trip absent, zero, and positive branch lengths as whole MAT/graph values.
- Verify that downstream algorithms retain their existing finite, non-negative validation where branch lengths are interpreted as evolutionary distances.
- Run the full lint and test suite.

## Related issues

Source: [kb/issues/M-io-usher-missing-branch-length-becomes-zero.md](../issues/M-io-usher-missing-branch-length-becomes-zero.md)
