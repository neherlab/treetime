# Preserve UShER zero branch lengths

Retain the distinction between an explicit zero branch length and an absent UShER branch length.

## Required changes

- Decode presence without filtering numeric zero.
- Encode `Some(0.0)` as an explicit zero and `None` as absence.
- Apply the common finite/non-negative branch-length validation at the format boundary.

## Validation

- Round-trip absent, zero, and positive branch lengths as whole MAT/TreeIR values.
- Reject negative, NaN, and infinite inputs with field context.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-io-usher-zero-branch-length-collapsed.md](../issues/M-io-usher-zero-branch-length-collapsed.md)
