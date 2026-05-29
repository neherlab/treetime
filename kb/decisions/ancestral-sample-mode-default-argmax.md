# Ancestral SampleMode default: Argmax

v1 defaults `SampleMode` to `Argmax` (deterministic) rather than `Root` (stochastic sampling at root), which is what augur always passes to v0 via `sample_from_profile='root'`.

## v0 behavior

v0's `prof2seq()` in `packages/legacy/treetime/treetime/seq_utils.py` accepts `sample_from_prof` as a parameter. Augur hardcodes `sample_from_profile='root'` in `packages/legacy/treetime/treetime/treetime.py`, which samples from the posterior CDF at the root node and uses argmax everywhere else. The RNG is seeded from `--seed` on the `TreeAnc` instance.

## v1 behavior

v1 defines `SampleMode` in `packages/treetime/src/ancestral/sample.rs` with variants `Argmax`, `Root`, `All`. The default is `Argmax`. The CLI `--sample-from-profile` flag (pending wiring) will allow selecting `Root` to match augur behavior.

## Rationale

Deterministic default ensures reproducible output without requiring a seed parameter. Stochastic sampling only affects positions where the root posterior has near-equal probabilities for multiple states. Users who need augur-compatible behavior can opt in via `--sample-from-profile=root --seed=N`. The default matches v1's existing deterministic `argmax_first()` behavior.
