# Fitch polytomy recurrence does not minimize arbitrary-degree mutations

## Problem

The Fitch backward pass intersects every child state set at a multifurcation and takes their union when the common intersection is empty. This is exact for binary nodes. For three or more children, one union event does not encode how many child edges must change, so the retained state set can include assignments with different mutation costs.

For three fixed children `A`, `A`, and `C`, the current recurrence returns `{A,C}`. Assigning the parent `A` costs one mutation, while assigning it `C` costs two. The forward pass can therefore select a non-minimum assignment depending on the parent state above the polytomy. TreeTime v0 uses the same all-child intersection/union rule, so correcting this changes parity and reconstructed ancestral states.

## Impact

- Binary trees remain exact under the standard unordered-character Fitch recurrence.
- Multifurcating trees can receive non-minimal ancestral assignments and mutation placement.
- Sparse marginal initialization inherits the selected Fitch assignment.

## Evidence

- v1 recurrence: [`packages/treetime/src/ancestral/fitch_sub.rs`](../../packages/treetime/src/ancestral/fitch_sub.rs) (`resolve_informative_positions_backward`)
- v0 recurrence: [`packages/legacy/treetime/treetime/treeanc.py`](../../packages/legacy/treetime/treetime/treeanc.py) (`_fitch_ancestral`)
- Existing multifurcation characterization: [`packages/treetime/src/ancestral/__tests__/test_fitch.rs`](../../packages/treetime/src/ancestral/__tests__/test_fitch.rs) (`test_fitch_polytomy`)

## Required fix

Use an exact finite-state cost recurrence for arbitrary node degree, retain every minimum-cost parent state, and validate the result against exhaustive assignments on generated multifurcating trees. This requires an approved scientific-output and v0-parity change.
