# Standardize mutations and tree-format projection

Make zero-based positions and typed mutation variants the single application contract, then perform checked conversion only at external format boundaries.

## Required changes

1. Define one mutation representation covering substitutions, insertions, and deletions with zero-based positions.
2. Add canonical mutation-string methods to that type. Remove duplicated mutation parsers and string assembly from TreeIR and command output code.
3. Convert one-based formats using checked addition/subtraction and checked integer narrowing. Reject zero, negative, and overflow values with contextual errors.
4. Project mutation variants only through existing format mappings that are verified as lossless.
5. Do not define UShER or PhyloXML loss policy in this ticket. Those mappings remain blocked by their related decision issues.

For an internal position $p_0$, serialize the external position as $p_1=p_0+1$. Parse only $p_1\ge1$ and recover $p_0=p_1-1$.

## Validation

- Parameterized unit tests for all mutation variants and boundary positions.
- Property tests for internal/string/format round trips.
- Boundary rejection cases for external zero and negative coordinates, values above `i32::MAX` for MAT, and reference lookups beyond sequence length.
- Independent fixtures for every verified lossless format mapping.
- Regression tests proving unresolved UShER and PhyloXML variants cannot be silently omitted by code touched here.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-core-mutation-representation-and-format-projection-inconsistent.md](../issues/M-core-mutation-representation-and-format-projection-inconsistent.md)
- Blocked format policies: [kb/issues/M-io-usher-mat-mutation-loss-is-implicit.md](../issues/M-io-usher-mat-mutation-loss-is-implicit.md)
- Blocked format policies: [kb/issues/N-io-phyloxml-mutation-property-contract-undecided.md](../issues/N-io-phyloxml-mutation-property-contract-undecided.md)
