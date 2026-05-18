# Sparse variable-site alphabet mismatch

v1 classifies variable sites using [`alphabet.is_ambiguous()`](../../packages/treetime/src/partition/sparse.rs#L43) while v0 uses `alphabet_gapN` (the alphabet extended with gap `-` and unknown `N` characters) in `treeanc.py`. The two approaches handle ambiguous characters differently, causing positions to be classified differently between v0 and v1. This affects which positions receive full marginal probability vectors vs shortcut computation.

Contributes to [dense-sparse log-likelihood divergence](ancestral-dense-sparse-log-likelihood-divergence.md)
and may affect v0/v1 numerical parity for alignments with gaps or ambiguous characters.

## Related issues

- Source: [M-ancestral-sparse-alphabet-mismatch.md](../issues/M-ancestral-sparse-alphabet-mismatch.md) -- delete after full resolution
- [M-ancestral-dense-sparse-divergence.md](../issues/M-ancestral-dense-sparse-divergence.md) -- dense-sparse log-likelihood divergence (related root cause)
