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

The TBP/calendar coordinate mismatch (formerly H-timetree-coalescent-coordinate-mismatch) has been fixed. Coalescent contributions now operate in calendar time and have a real effect on node times, making this skipped-initial-Tc issue observable.
