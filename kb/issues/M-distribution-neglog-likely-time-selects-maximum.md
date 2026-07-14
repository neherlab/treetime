# Neg-log likely time selects the least-likely ordinate

`Distribution<Y>::likely_time()` delegates both `Plain` and `NegLog` functions to helpers that select the maximum stored ordinate:

- [`packages/treetime-distribution/src/distribution_core/distribution.rs#L67-L83`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L67-L83)
- [`packages/treetime-distribution/src/distribution_core/function.rs#L248-L259`](../../packages/treetime-distribution/src/distribution_core/function.rs#L248-L259)
- [`packages/treetime-distribution/src/distribution_core/formula.rs#L66-L78`](../../packages/treetime-distribution/src/distribution_core/formula.rs#L66-L78)

For `DistributionPlain`, the maximum ordinate is the likelihood peak. For `DistributionNegLog`, defined at [`packages/treetime-distribution/src/distribution_core/distribution.rs#L628`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L628), lower ordinates represent higher likelihoods. Selecting the maximum therefore returns the least-likely sampled time for `Function` and `Formula` variants.

## Impact

The public generic API returns the wrong result for a valid `DistributionNegLog`. Current main inference paths commonly normalize or convert distributions before selecting a likely time, which limits observed exposure, but the type itself makes no such precondition.

## Required behavior

- Select the time whose ordinate represents the maximum plain-space likelihood under the active `YAxisPolicy`.
- Preserve existing `Empty`, `Point`, and `Range` behavior.
- Test matching `Plain` and `NegLog` functions and formulas so that a plain-space peak and its neg-log trough identify the same time.

## Related issues

- [M-timetree-coalescent-missing-leaf-and-root-contributions.md](M-timetree-coalescent-missing-leaf-and-root-contributions.md) — coalescent objective terms use `DistributionNegLog`, so consumers must respect its ordinate semantics.

## Related tickets

- [kb/tickets/distribution-make-likely-time-policy-aware.md](../tickets/distribution-make-likely-time-policy-aware.md)
