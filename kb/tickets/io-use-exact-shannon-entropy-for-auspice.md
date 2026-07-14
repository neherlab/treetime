# Use exact Shannon entropy for Auspice

Replace the perturbed hand-written entropy formula with the maintained ndarray-stats implementation.

## Required changes

- Consume one borrowed profile row.
- Use `EntropyExt::entropy()` without adding an epsilon to probabilities.
- Propagate empty, non-finite, and invalid-profile errors with node/trait context.

## Validation

- Deterministic, uniform, zero-containing, empty, and invalid distributions.
- Independent analytical expected values.
- Auspice whole-document projection test.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-io-auspice-entropy-perturbs-shannon-definition.md](../issues/M-io-auspice-entropy-perturbs-shannon-definition.md)
