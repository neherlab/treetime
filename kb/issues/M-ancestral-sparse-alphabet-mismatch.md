# Sparse variable-site alphabet mismatch

v1 classifies variable sites using [`alphabet.is_ambiguous()`](../../packages/treetime/src/partition/payload/sparse.rs#L43) while v0 uses `alphabet_gapN` (the alphabet extended with gap `-` and unknown `N` characters) in `treeanc.py`. The two approaches handle ambiguous characters differently, causing positions to be classified differently between v0 and v1. This affects which positions receive full marginal probability vectors vs shortcut computation.

Contributes to [dense-sparse log-likelihood divergence](M-ancestral-dense-sparse-divergence.md) and may affect v0/v1 numerical parity for alignments with gaps or ambiguous characters.
