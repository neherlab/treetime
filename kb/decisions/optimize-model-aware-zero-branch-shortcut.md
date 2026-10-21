# Model-aware zero-branch shortcut

## Deviation

v1 restricts the `is_zero_branch_optimal()` derivative-at-zero shortcut to GTR model families with proven unimodal branch-length likelihood (JC69, F81, binary models). For models with multiple distinct nonzero eigenvalues (K80, HKY85, T92, TN93, JTT92, inferred GTR), the shortcut is bypassed and Newton/grid search evaluates the full likelihood surface.

v0 does not have this shortcut at all -- it uses Brent's method (`scipy.optimize.minimize_scalar`) for all branch length optimization regardless of model.

## Rationale

Dinh & Matsen (2017) prove that for models with one distinct nonzero eigenvalue (JC69, F81, binary), the one-dimensional phylogenetic likelihood $L(t)$ has at most one stationary point on $(0, \infty)$. A negative derivative at $t = 0$ guarantees zero is the global maximum.

For K80 (Kimura 2-parameter), the same paper constructs a counterexample with two local maxima and a local minimum. A negative derivative at $t = 0$ only proves zero is a local maximum -- a better positive-length maximum may exist. The same multimodality applies to HKY85, TN93, and general GTR where eigenvalue ratios differ.

The unrestricted shortcut could collapse a branch to zero even when a positive branch length has higher likelihood.

## Implementation

- `fn GTR::certified_decay_rate()` [packages/treetime/src/gtr/gtr.rs#L130-L132](../../packages/treetime/src/gtr/gtr.rs#L130-L132): the decay rate shared by every nonzero eigenvalue, `None` for models where $L(t)$ is not proven unimodal
- `fn GTR::new()` sets the unit decay rate only when all exchangeabilities are equal and positive [packages/treetime/src/gtr/gtr.rs#L100](../../packages/treetime/src/gtr/gtr.rs#L100) (`fn has_equal_positive_exchangeabilities()` [packages/treetime/src/gtr/gtr.rs#L256](../../packages/treetime/src/gtr/gtr.rs#L256)), which holds for the F81 family: JC69, F81, binary models, TN93 selected by name (default `kappa1 = kappa2 = 1`), and custom, inferred, or mugration models with equal exchangeabilities. K80, HKY85, and T92 selected by name use the default kappa 0.1 and are not certified
- `fn is_zero_branch_optimal()` [packages/treetime/src/optimize/zero_boundary.rs#L21-L40](../../packages/treetime/src/optimize/zero_boundary.rs#L21-L40): applies the derivative shortcut only when `fn shared_certified_decay_rate()` [packages/treetime/src/optimize/zero_boundary.rs#L120-L127](../../packages/treetime/src/optimize/zero_boundary.rs#L120-L127) finds the same certified decay rate in every contribution

## Impact

- JC69 and F81 (default and most common models): no change
- K80/HKY85/T92 (default kappa 0.1), JTT92, and inferred or custom models with unequal exchangeabilities: shortcut bypassed, branches that would have been set to zero now go through Newton/grid search
- TN93 selected by name (default kappa 1) and any other model with equal exchangeabilities: certified, shortcut applies
- Minor performance cost for trees with many near-zero edges under non-simple models

## References

- Dinh V, Matsen FA IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Ann. Appl. Probab._ 27(3):1646-1677. https://doi.org/10.1214/16-AAP1240
