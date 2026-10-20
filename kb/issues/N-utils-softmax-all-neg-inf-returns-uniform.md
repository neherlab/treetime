# `softmax_with_log_norm_inplace` returns a uniform distribution when every state is impossible

## Summary

`softmax_with_log_norm_inplace` in `packages/treetime-utils/src/array/softmax_with_log_norm.rs` returns a uniform distribution with `log_norm = -inf` when every input entry is $-\infty$. The unit test `test_softmax_with_log_norm_degenerate` in `packages/treetime-utils/src/array/__tests__/test_softmax_with_log_norm.rs` fixes this behavior in place without stating a reason.

## Consequence

Marginal reconstruction normalizes each site's log-likelihood row with this function in `packages/treetime/src/partition/marginal/shared/normalize.rs` and in `packages/treetime/src/partition/marginal/sparse/message.rs`. A site where every state has zero likelihood therefore receives a uniform profile, and the summed log-likelihood becomes $-\infty$. The uniform profile does not show that the site is impossible under the model, so downstream consumers of the profile cannot distinguish it from a site with no information.

## v0 behavior

v0 `normalize_profile()` with `log=True` subtracts the row maximum and exponentiates [packages/legacy/treetime/treetime/seq_utils.py#L296-L306](../../packages/legacy/treetime/treetime/seq_utils.py#L296-L306). For an all-$-\infty$ row it computes $\exp(-\infty - (-\infty)) = \exp(\text{NaN})$, so the profile and the log offset are NaN. v1 differs from v0 on this input: it returns a uniform profile from the branch at [packages/treetime-utils/src/array/softmax_with_log_norm.rs#L25-L28](../../packages/treetime-utils/src/array/softmax_with_log_norm.rs#L25-L28).

> [!IMPORTANT]
> **Decision required.** The contract for an all-impossible row is undecided. v0 gives NaN, v1 gives a uniform profile, and neither names the site. Options:
>
> - **Keep the uniform fallback**: record it as intended behavior and rely on the $-\infty$ log-likelihood to signal the condition
> - **Return an error**: report the site as impossible under the model, so the caller can name the site and the partition
> - **Return NaN probabilities**: match v0 and propagate the undefined distribution, as the function already does for NaN input
