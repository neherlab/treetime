# --coalescent-opt alone skips initial Tc pass

When only `--coalescent-opt` is specified (without a fixed `--coalescent=<N>`),
the first iteration runs without the coalescent prior. The coalescent model is
applied only in subsequent iterations after the initial Tc is optimized.

v0 always enables the coalescent from the first iteration when `--coalescent-opt`
is used, starting with an initial Tc estimate.

## Impact

The first iteration produces node times without coalescent regularization. The
Tc optimization then operates on these unregularized times, producing a
different Tc estimate than v0. Subsequent iterations use the coalescent prior but
start from a different initial state.

## Related issues

Coalescent contributions operate in calendar time, so skipping the initial Tc pass changes the node times that seed later iterations.

This issue compounds with other coalescent-prior gaps. Any missing contribution in the backward pass or convergence metric is magnified when the first iteration already starts from an unregularized state.
