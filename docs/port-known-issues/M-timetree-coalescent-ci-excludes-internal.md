# Coalescent and skyline CI excludes internal nodes

When `--coalescent-skyline --time-marginal=only-final` (or any fixed Tc
coalescent with `--time-marginal=only-final`) are combined, the
`confidence_intervals.tsv` file contains only tips. Internal nodes are excluded
because their time distributions are `None` before marginal reconstruction runs.

## Example

flu/h3n2/20 with `--coalescent-skyline --time-marginal=only-final`: 19 entries
(tips only, no internal nodes).

## Root cause

Compound of two issues:

1. Internal node dates missing at scale (now fixed) dropped internal node times
   to `None`
2. The marginal reconstruction at `only-final` cannot recover these nodes because
   the time distributions are already `None`

CI estimation is impossible for internal nodes when using coalescent (fixed Tc)
or skyline models.

## Related issues

- Internal node dates missing at scale (now fixed): was the primary cause
- With the scale fix applied, this issue may be partially resolved. Needs revalidation.
