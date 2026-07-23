# Timetree likelihood API hides the log domain

The convergence API names functions `compute_*_likelihood` and serialized fields `lh_seq`, `lh_pos`, `lh_coal`, and `lh_total`, although every value is a log-likelihood and totals are formed by addition [`packages/treetime/src/timetree/convergence/likelihood.rs#L12`](../../packages/treetime/src/timetree/convergence/likelihood.rs#L12) [`packages/treetime/src/timetree/convergence/metrics.rs#L10`](../../packages/treetime/src/timetree/convergence/metrics.rs#L10).

The broader codebase already distinguishes `log_lh` and `log_likelihood`, including `graph_log_lh()` and Poisson indel log-likelihoods. The convergence surface erases that distinction at a public and serialized boundary, making aggregation, comparison, and display semantics easy to misread.

## Observable mismatch

- Sequence likelihood delegates to `graph_log_lh()`.
- Positional likelihood sums logarithms of probabilities.
- Coalescent likelihood delegates to a total log-likelihood operation.
- `lh_total` adds components, which is valid in log space and would be invalid for ordinary likelihoods.

The names therefore permit a caller to exponentiate, average, compare, or display the values under the wrong mathematical interpretation. The serialized field names propagate that ambiguity beyond one module.

## Required vocabulary

Rename functions, fields, documentation, and output labels to state `log_likelihood` or the established `log_lh` abbreviation consistently. v1 has not shipped, so compatibility aliases would preserve ambiguity without serving an existing consumer.

This rename does not define a new convergence criterion or change which likelihood components are present; those behaviors are tracked separately.

## Validation

- Unit tests verify each component against its current log-domain oracle.
- Total composition remains the sum of available log-likelihood components.
- Serde and validation-output tests assert the renamed schema explicitly.
- Repository search finds no convergence API that labels a log-likelihood as an ordinary likelihood.

## Related issues

- [M-timetree-convergence-metric-deficiencies.md](M-timetree-convergence-metric-deficiencies.md)
- [N-timetree-convergence-metric-excludes-coalescent.md](N-timetree-convergence-metric-excludes-coalescent.md)
