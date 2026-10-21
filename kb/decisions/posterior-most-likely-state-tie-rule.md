# The most likely state is the first state within 1e-10 of the maximum

Every decision of a most likely state over a posterior takes the first state, in the fixed state order, whose probability is at most $10^{-10}$ below the largest probability of the posterior. v0 takes the first exact maximum, so v1 diverges from v0 where two states of a posterior differ by $10^{-10}$ or less.

**Type**: Scientific output; divergence from v0.

**Status**: Approved on 2026-10-21.

## Rule

- `most_likely_state()` in [packages/treetime/src/partition/marginal/posterior_states.rs](../../packages/treetime/src/partition/marginal/posterior_states.rs) applies the rule with `POSTERIOR_TIE_TOLERANCE` of [packages/treetime/src/constants.rs](../../packages/treetime/src/constants.rs), through `argmax_first_within()` of `treetime-utils`
- Every most likely state over a posterior comes from this function: the state of a trait (`StatePosterior::new()`), the reconstructed sequences of the dense and the sparse representation, the states of branch mutations, the imputed states of leaves, and the root states that the GTR fit counts from the marginal posteriors
- Two probabilities are tied when they differ by at most $10^{-10}$ in absolute value. All posteriors are probabilities on one scale from 0 to 1. A tolerance relative to the maximum was rejected
- The fixed state order is the order of the canonical states of the alphabet (`A`, `C`, `G`, `T` for nucleotides) and the sorted state names for traits (`StateAlphabet`)
- The sparse imputation of a leaf normalizes the posterior restricted to the observed character before it applies the rule, so the dense and the sparse representation compare the same probabilities

## Reference behavior

`prof2seq()` of v0 takes `numpy.argmax` of the normalized profile (`packages/legacy/treetime/treetime/seq_utils.py` line 271), which returns the first state of the exact maximum. The ancestral and the mugration reconstruction of v0 use it.

## Why

States that are equal in exact arithmetic differ in the last bits of their computed probabilities, and these bits depend on the BLAS, LAPACK, and math library of the platform. Under the exact maximum, the reconstructed state, the mutations, the branch labels, and the Auspice colorings of such a node change from platform to platform. Example: `mugration` on `data/dengue/100` with the attribute `country` gave node `NODE_0000040` the probabilities 0.3916531565060838 (`Samoa`) and 0.39165315650588367 (`American Samoa`) with one OpenBLAS build, and 0.39165315650834764 (`American Samoa`) and 0.3916531565081366 (`Samoa`) with another (commits `2e275920` and `4f116c6d`). Under the rule, both builds give `American Samoa`, the first of the two names, and `--ambiguity-threshold=0` reports `American Samoa|Samoa` for this node only. The reproducibility target is agreement within tolerances ([kb/decisions/reproducibility-across-platforms-within-tolerances.md](reproducibility-across-platforms-within-tolerances.md)), and a decision over computed values needs this rule to meet it.

## Limits

The rule removes the dependence on rounding only for differences up to $10^{-10}$. Two states whose probabilities differ by more than that on one platform and by less on another still give different states. The optional ambiguity report names all states near the most likely state ([kb/decisions/ambiguity-report-of-plausible-states.md](ambiguity-report-of-plausible-states.md)).
