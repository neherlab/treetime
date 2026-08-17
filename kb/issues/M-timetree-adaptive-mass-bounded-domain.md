# Branch grid uses ad-hoc peak multiples instead of a mass-bounded domain

`create_simple_grid()` sets the grid extent from fixed multiples of the peak branch length and one mutation, not from a target probability mass. The proposal replaces this with a mass-bounded domain, re-derived once per operation, with the boundary laws carrying the remaining tail mass.

## Mechanism

The current extent is `min_bl = one_mutation * 0.01` to `max_bl = min(max(center * 5, one_mutation * 10), 5.0)` [`packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L78-L97`](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L78-L97). The `5 * center` multiple is the tight case: for a one-mutation branch it leaves about four percent of the Gamma mass outside the window. No quantile-based re-windowing exists; the only quantile code computes confidence intervals, not grid domains.

## Required behavior

Derive the domain from mass, with both boundary regions contributing in closed form:

```text
domain = [ max(hard_lo, q(eps)),  min(hard_hi, q(1 - eps)) ]
```

- Suppress the eps-trim on any side where the mode sits on the hard bound, so a mode on a hard edge is never trimmed away.
- Use a single named epsilon per side (proposed 5e-4), not a tuned multiple.
- Apply the same rule at construction (replacing the ad-hoc extent) and after every operation.
- Use a fixed point count with adaptive spacing and a resolution floor, so multiplying a narrow distribution by a wide one cannot ratchet the spacing coarser than the narrow operand had.

The soft-tail and hard-approach mass integrals are closed form, so the quantile integral is cheap.

## Impact

The finite window was never the correctness bug (the boundary semantics were), but a mass-bounded domain removes the tuned `5 *` multiple and makes the grid a deliberate mass target rather than a support. This is proposal step 7. Its prerequisites are in place: the common-grid backward pass is implemented, and the messages now carry a fitted integrable `Linear` soft tail whose mass integral is closed form ([kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md)).

## Related

- [M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md): the grid-spacing and golden-master acceptance concern for the same grid.
- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): Part D.
