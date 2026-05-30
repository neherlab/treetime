# Ancestral SampleMode default: Argmax

v1 defaults `SampleMode` to `Argmax` (deterministic) rather than `Root` (stochastic sampling at root), which is what augur always passes to v0 via `sample_from_profile='root'`.

## v0 behavior

v0's `prof2seq()` in `packages/legacy/treetime/treetime/seq_utils.py` accepts `sample_from_prof` as a parameter. Augur hardcodes `sample_from_profile='root'` in `packages/legacy/treetime/treetime/treetime.py`, which samples from the posterior CDF at the root node and uses argmax everywhere else. The RNG is seeded from `--seed` on the `TreeAnc` instance.

## v1 behavior

v1 defines `SampleMode` in `packages/treetime/src/ancestral/sample.rs` with variants `Argmax`, `Root`, `All`. The default is `Argmax`. The `ancestral` command exposes `--sample-from-profile=<argmax|root|all>`, threading the chosen mode and a `--seed`-seeded RNG through `ancestral_reconstruction_marginal` into the sparse and dense reconstruction paths. `Root` samples the posterior at the root only (argmax elsewhere), matching augur's `sample_from_profile='root'`; `All` samples at every node. Sampling applies only to marginal reconstruction; the forward-pass and convergence reads stay deterministic regardless of mode.

## Rationale

Deterministic default ensures reproducible output without requiring a seed parameter. Stochastic sampling only affects positions where the posterior has near-equal probabilities for multiple states. Users who need augur-compatible behavior can opt in via `--sample-from-profile=root --seed=N`. The default matches v1's existing deterministic `argmax_first()` behavior.

`Root`/`All` are reproducible within v1 under a fixed `--seed` (single-threaded preorder draw order, seeded `StdRng`), but are not byte-identical to v0: v0 draws from numpy's PCG64 while v1 uses `StdRng`. Only `Argmax` reproduces v0 output exactly, which is why the golden-master and Python-parity tests run in argmax mode.
