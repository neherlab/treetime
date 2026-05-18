# Coalescent leaf survival formula sign convention requires investigation

## Summary

The coalescent contribution functions use different sign conventions for leaf vs internal node survival integrals, and the return type name (`DistributionNegLog`) contradicts the actual value convention (log-probability). Requires verification against v0 and theory.

## Details

`packages/treetime/src/coalescent/contributions.rs`

The docstring states:

- Leaves: exp(I(t)) (survival probability, positive exponent)
- Internal nodes: exp(-I(t)) (negative exponent)

Line 55 states: "The sign convention follows Python v0 which stores -I(t) in neg-log space."

The return type is `DistributionNegLog`, which by naming convention implies the stored values are negated log-probabilities (positive values = low probability). The actual values follow log-probability convention (more negative = less likely), contradicting the type name.

## Open questions

1. Is the sign difference between leaf and internal contributions intentional or a bug?
2. Does the `DistributionNegLog` type store neg-log values or log values? The name and usage disagree.
3. Does v0 use the same sign convention, or was the sign inherited from a misreading of v0?
4. What is the theoretical justification for different signs on leaf vs internal survival?

## Required investigation

- Trace the v0 Python implementation to verify sign conventions for leaf and internal coalescent contributions
- Verify against Volz 2012 or Sagulenko 2018 coalescent formulation
- Determine whether `DistributionNegLog` is misnamed or misused
- If the sign is wrong: quantify impact on coalescent prior and time estimates

## Impact

- If sign is wrong, coalescent prior contribution is inverted for leaves or internals
- Incorrect coalescent prior biases time estimates toward or away from the present
- The type-name confusion makes the code harder to audit for correctness
