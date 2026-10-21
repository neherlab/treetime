# The reconstructed state of mugration at a tied node follows the rounding of the platform

When two states of a node have the same marginal posterior, the reconstructed state of `mugration` is the one whose probability the last bits of the linear algebra make larger. The same inputs then give different states, branch labels, and Auspice colorings on platforms whose BLAS or LAPACK round differently.

## Evidence

- `dev/smoke` of commit `2e275920` against `4f116c6d` (OpenBLAS develop `de5d4914` against the OpenBLAS 0.3.34 release, `Haswell` kernels on both): on `mugration/dengue/100/country`, node `NODE_0000040` has the posteriors `Samoa` 0.3916531565060838 and `American Samoa` 0.39165315650588367 in the baseline, and `American Samoa` 0.39165315650834764 and `Samoa` 0.3916531565081366 in the current build. The two values differ by 5e-13 relative, every other value of the run by at most 4.6e-9, and the reconstructed state of the node changes from `Samoa` to `American Samoa`, with the branch labels of `NODE_0000040`, `HM582105`, and `OK469349`
- `StatePosterior::new()` [packages/treetime/src/partition/marginal/discrete/state_posterior.rs#L20-L22](../../packages/treetime/src/partition/marginal/discrete/state_posterior.rs#L20-L22) takes the first maximum of the posterior profile (`argmax_first`), which is deterministic only for values that are equal bit for bit. Every output reads the state from it

> [!IMPORTANT]
> **Decision required.** The rule for states whose posteriors agree within rounding:
>
> - take the first state in alphabet order among the states within a relative tolerance of the maximum, so the choice is the same on every platform
> - keep the exact maximum and accept platform-dependent states at ties

## Validation

- `mugration` on `data/dengue/100` gives the same reconstructed state at `NODE_0000040` with the OpenBLAS kernels of different CPUs and in web2
