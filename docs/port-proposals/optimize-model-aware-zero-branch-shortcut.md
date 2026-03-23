# Model-aware zero-branch shortcut

Restrict the `is_zero_branch_optimal()` shortcut to model families with a proven uniqueness result for the one-dimensional branch-length likelihood. For models where the likelihood can have multiple local optima, skip the shortcut and fall through to Newton/grid search.

## Scientific motivation

Dinh and Matsen (2017) prove that for JC69, F81, and all binary models, the one-dimensional phylogenetic likelihood function `L(t)` has at most one stationary point on `(0, ∞)`. For these models, a negative derivative at `t = 0` guarantees zero is the global maximum on `[0, ∞)`.

For K80 (Kimura 2-parameter), the same paper demonstrates that `L(t)` can have two local maxima and a local minimum. A negative derivative at `t = 0` proves only that zero is a local maximum, not that no better positive-length maximum exists. The same multimodality applies to HKY85, TN93, and general GTR models where eigenvalue ratios differ from JC69/F81.

## Current behavior

`is_zero_branch_optimal()` applies the derivative-sign shortcut to all models. For K80/HKY85/TN93/GTR, this can collapse a branch to zero even when a positive branch length has higher likelihood.

## Proposed change

Check the active GTR model before applying the shortcut. If the model belongs to a family with proven unimodality (JC69, F81, binary models), use the shortcut. For all other models, return false and let Newton/grid search evaluate the full likelihood surface.

## Expected impact

For JC69 and F81 (the default and most common models): no change.
For K80/HKY85/TN93/GTR: the shortcut is bypassed. Branches that the shortcut would have set to zero now go through Newton/grid search, which correctly finds the global maximum. Minor performance cost for trees with many near-zero edges under non-simple models.

## Validation plan

- Construct a K80 model with coefficients producing two local maxima in `L(t)` where `L'(0) < 0` but the global maximum is at `t > 0`
- Verify the shortcut returns true (current behavior, incorrect) and Newton/grid finds the correct positive maximum
- After the fix, verify the shortcut returns false for K80 and the optimizer finds the correct maximum

## References

- Dinh V, Matsen FA IV. "The shape of the one-dimensional phylogenetic likelihood function." Annals of Applied Probability 27(3), 2017. doi:10.1214/16-AAP1240. PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC10153603/
