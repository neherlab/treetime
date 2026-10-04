# `timetree` takes the most likely root sequence where v0 samples it

v0 `timetree` samples the root sequence from its posterior in every ancestral reconstruction of the run. v1 `timetree` always takes the most likely state at the root. No decision records this difference, so it is an unapproved v0 divergence until the team decides.

## What v0 does

`TreeTime._run()` passes `sample_from_profile='root'` to every reconstruction of the run: the first one, the one after rerooting, the ones after polytomy resolution, and the one in every iteration ([packages/legacy/treetime/treetime/treetime.py#L211-L217](../../packages/legacy/treetime/treetime/treetime.py#L211-L217), [packages/legacy/treetime/treetime/treetime.py#L345-L348](../../packages/legacy/treetime/treetime/treetime.py#L345-L348)). The value is hardcoded and no command-line flag changes it. Both reconstruction methods sample at the root and take the most likely state at all other nodes:

- Marginal: `_ml_anc_marginal()` ([packages/legacy/treetime/treetime/treeanc.py#L785-L792](../../packages/legacy/treetime/treetime/treeanc.py#L785-L792))
- Joint: `_ml_anc_joint()` ([packages/legacy/treetime/treetime/treeanc.py#L1023-L1033](../../packages/legacy/treetime/treetime/treeanc.py#L1023-L1033))

`prof2seq()` draws one uniform number per alignment column from the generator seeded by `--rng-seed` and picks the state by inverse CDF ([packages/legacy/treetime/treetime/seq_utils.py#L266-L271](../../packages/legacy/treetime/treetime/seq_utils.py#L266-L271)). Without `--rng-seed`, the seed comes from entropy and is not reported. Augur `refine` runs through the same `TreeTime.run()`, so its node data has a sampled root sequence too.

## What v1 does

v1 `timetree` writes the sequences and mutations of the marginal reconstruction with the most likely state at every node ([packages/treetime/src/timetree/pipeline.rs#L549-L555](../../packages/treetime/src/timetree/pipeline.rs#L549-L555), `fn node_sequence()` of `enum MarginalReconstruction` in `packages/treetime/src/partition/marginal/reconstruction.rs`). Sampling from the posterior exists only in the `ancestral` command (`--sample-from-profile`). [kb/decisions/ancestral-sample-mode-default-argmax.md](../decisions/ancestral-sample-mode-default-argmax.md) approves the most-likely-state default for `ancestral` only.

## Impact

- At columns where the root posterior has no dominant state, v1 and v0 can write a different root state, and therefore different mutations on the branches below the root, in the reconstructed sequences, augur node-data JSON, Auspice JSON, and UShER MAT outputs
- v0 output at those columns changes with the seed; v1 output is deterministic
- In v0 `joint` branch-length mode, the sampled root also changes the branch lengths below the root, and with them the node dates

## Potential solutions

- Approve the most likely root state for `timetree` as an intentional change, and extend the decision to cover it
- Match v0: sample the root in every `timetree` reconstruction. This adds a random step to every `timetree` run, so it must follow [kb/decisions/cli-seed-on-commands-with-random-steps.md](../decisions/cli-seed-on-commands-with-random-steps.md): draw from `get_random_number_generator()` with the run seed that `SeedArgs::resolve()` returns (the `--seed` value, or a drawn seed that the run logs). `SeedArgs::resolve()` must then name root sampling as a random step whenever an alignment is present, not only polytomy resolution
- Expose the root sampling as an option of `timetree`, with the most likely state as the default

## Related issues

- [M-clock-alignment-ignored.md](M-clock-alignment-ignored.md): v0 `clock` reaches the same root sampling through `TreeTime.run()` when it uses the alignment
