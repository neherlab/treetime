# The trait entropy is the exact Shannon entropy

The entropy that the trait outputs write for a node is the Shannon entropy of its state posterior, $H = -\sum_i p_i \ln p_i$ with $0 \ln 0 = 0$. Augur adds $10^{-12}$ inside the logarithm, so v1 diverges from augur by at most $n \cdot 10^{-12}$ for $n$ states.

**Type**: Scientific output; divergence from augur.

**Status**: Approved on 2026-10-21.

## Rule

`StatePosterior::new()` ([packages/treetime/src/partition/marginal/discrete/state_posterior.rs](../../packages/treetime/src/partition/marginal/discrete/state_posterior.rs)) computes the entropy with `EntropyExt::entropy()` of ndarray-stats, which counts a zero probability as zero. The augur node data and the Auspice JSON read it from there ([kb/decisions/trait-probabilities-written-at-full-precision.md](trait-probabilities-written-at-full-precision.md)).

## Reference behavior

`augur traits` computes `S = -np.sum(pdis*np.log(pdis+TINY))` with `TINY = 1e-12` (`augur/traits.py` lines 13 and 87 at revision `58f329042e7fc631ca275cf2df44c8720f80949f`). The offset came with the first augur mugration code (commit `a57201a6`, 2018-04-29) without a stated reason. It keeps numpy from computing `0 * log(0)` as NaN for the exact zeros of a tip with a known state. v0 computes no trait entropy.

## Why

The offset changes every term: the entropy of each node is lower by about $k \cdot 10^{-12}$ for $k$ states with a nonzero probability, and every tip with a known state gets $-1.00009 \cdot 10^{-12}$, a negative entropy, which is impossible. The convention $0 \ln 0 = 0$ also avoids NaN, without changing any term, and the outputs write the entropy at full precision, where the offset is visible.
