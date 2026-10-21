# The most likely state at a tie follows the rounding of the platform

When two states of a posterior are equal in exact arithmetic, the reconstructed state is the one that the last bits of the linear algebra make larger. The same inputs then give different reconstructed states, mutations, branch labels, and Auspice colorings on platforms whose BLAS, LAPACK, or math library round differently. Exact ties appear where the data and the model are symmetric in two states.

## Evidence

- `dev/smoke` of commit `2e275920` against `4f116c6d` (OpenBLAS develop `de5d4914` against the OpenBLAS 0.3.34 release, `Haswell` kernels on both): on `mugration/dengue/100/country`, node `NODE_0000040` has the posteriors `Samoa` 0.3916531565060838 and `American Samoa` 0.39165315650588367 in the baseline, and `American Samoa` 0.39165315650834764 and `Samoa` 0.3916531565081366 in the current build. The two values differ by 5e-13 relative, every other value of the run by at most 4.6e-9, and the reconstructed state of the node changes from `Samoa` to `American Samoa`, with the branch labels of `NODE_0000040`, `HM582105`, and `OK469349`
- Every decision of a most likely state takes the first maximum with `argmax_first()`, which is deterministic only for values that are equal bit for bit:
  - trait states: `StatePosterior::new()` [packages/treetime/src/partition/marginal/discrete/state_posterior.rs#L20-L21](../../packages/treetime/src/partition/marginal/discrete/state_posterior.rs#L20-L21)
  - dense sequences: `Resolve::Argmax` [packages/treetime/src/partition/marginal/sample.rs#L30](../../packages/treetime/src/partition/marginal/sample.rs#L30), and the states of branch mutations [packages/treetime/src/partition/marginal/dense/partition.rs#L95-L96](../../packages/treetime/src/partition/marginal/dense/partition.rs#L95-L96) and [packages/treetime/src/partition/marginal/dense/partition.rs#L214](../../packages/treetime/src/partition/marginal/dense/partition.rs#L214)
  - sparse sequences: [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L124](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L124) and [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L129](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L129)
  - the root states of the GTR fit: [packages/treetime/src/partition/marginal/shared/data.rs#L47](../../packages/treetime/src/partition/marginal/shared/data.rs#L47)

The reproducibility target is agreement within tolerances ([kb/decisions/reproducibility-across-platforms-within-tolerances.md](../decisions/reproducibility-across-platforms-within-tolerances.md)), so a decision over computed values needs a rule that gives the same state on every platform when the values are equal up to rounding.

## Decided

On 2026-10-21: one shared rule for every decision of a most likely state over posteriors; the ambiguity of tied states is reported as an option, in every output whose format can hold it, without limiting the feature to what the current formats support.

> [!IMPORTANT]
> **Decision required.** The criterion for "equal up to rounding". A tolerance relative to the maximum is rejected. Candidates: an absolute difference of probabilities, which share one scale from 0 to 1; a number of units in the last place of the larger value; or a criterion from the error bound of the computation.

> [!IMPORTANT]
> **Decision required.** The ambiguity report: which states it names (the ties of the rule above, or a statistical set such as the states within 0.5/n of the maximum, as RAxML-NG writes), and where: IUPAC codes in the sequence outputs, the set of states in the traits table, and extra fields in augur node data and the Auspice JSON.

## Validation

- `mugration` on `data/dengue/100` gives the same reconstructed state at `NODE_0000040` with the OpenBLAS kernels of different CPUs and in web2
- A symmetric two-state tree gives the same state for every permutation of the leaf order and on every kernel set
