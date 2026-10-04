# Seed handling: one seed per run, only on commands with a random step

v1 accepts `--seed` (alias `--rng-seed`) only on the commands that draw random numbers, `ancestral` and `timetree`. A run resolves one seed in the command layer and passes it to every random step. Without `--seed`, the command draws a seed and logs it whenever the run has a random step, so any run can be reproduced from its log.

## What v0 does

v0 adds `--rng-seed` to every subcommand parser in `packages/legacy/treetime/treetime/argument_parser.py` (for example `mugration` at [L402](../../packages/legacy/treetime/treetime/argument_parser.py#L402) and `clock` at [L443](../../packages/legacy/treetime/treetime/argument_parser.py#L443)) and passes the value to the NumPy generator of each `TreeAnc` instance. Without a seed, the generator is seeded from entropy and the seed is not reported. v0 draws random numbers for polytomy resolution (`packages/legacy/treetime/treetime/treetime.py`), profile sampling (`prof2seq()` in `packages/legacy/treetime/treetime/seq_utils.py`), sequence simulation, random GTR models, and the Fitch root tie-break (`packages/legacy/treetime/treetime/treeanc.py`).

## What v1 does

- The random steps of v1 are profile sampling in `ancestral` (`--sample-from-profile=root|all`) and polytomy resolution in `timetree` (`--resolve-polytomies`). The Fitch root tie-break is deterministic ([ancestral-fitch-deterministic-root-state.md](ancestral-fitch-deterministic-root-state.md)), so `clock` and `mugration` have no random step and do not accept `--seed`
- `struct SeedArgs` in `packages/app-commands/src/commands/shared/seed.rs` owns the flag. `SeedArgs::resolve()` returns the given seed or draws one, and logs `<step> is stochastic; seed <seed> (pass --seed to reproduce this run)` when the run has a random step
- The core takes a plain `u64` seed (`AncestralParams`, `AaParams`, `TimetreeParams`), and `get_random_number_generator()` in `packages/treetime-utils/src/sync/random.rs` takes a `u64`, so no core path reads entropy
- In `ancestral`, the nucleotide and amino-acid reconstructions of one run share the resolved seed. Each reconstruction starts its own generator from that seed, so `--seed=N` reproduces both

## Rule for every random step

Every operation that draws random numbers (sampling, random tie-breaking, stochastic search, simulation) follows the same contract, in core code and in commands:

- Draw only from a generator that `get_random_number_generator(seed)` returns, or from one passed down from it; there is no other source of randomness
- Take the seed as a required `u64` parameter of the core function or its params struct. The core never draws a seed and never reads entropy
- Resolve the seed once per run in the command, with `SeedArgs::resolve()`: the `--seed` value when given, otherwise a random seed. The command logs the seed whenever the run has a random step, so an unseeded run can be reproduced with `--seed=<logged value>`
- Pass that one seed to every random step of the run. A step may start its own generator from the run seed, as the nucleotide and amino-acid reconstructions of `ancestral` do, so that each step stays reproducible on its own
- A command that gains its first random step gains `--seed` by flattening `SeedArgs` into its arguments, and names the step in its `SeedArgs::resolve()` call so that the log line appears
- Tests seed every generator with a fixed value

Two v0 random steps are not ported yet and must follow this rule when they are: root sampling in `clock` with an alignment ([kb/issues/M-clock-alignment-ignored.md](../issues/M-clock-alignment-ignored.md)) and root sampling in `timetree` ([kb/issues/M-timetree-root-sequence-not-sampled.md](../issues/M-timetree-root-sequence-not-sampled.md)).

## Rationale

- A seed that no random step reads gives the user a false promise of reproducibility; [kb/issues/M-cli-flags-parsed-but-ignored.md](../issues/M-cli-flags-parsed-but-ignored.md) tracked the ignored `clock` and `mugration` flags
- Drawing one seed per run and logging it makes an unseeded stochastic run reproducible after the fact, which v0 does not offer
- v1 has no backward compatibility requirement, so a v0 invocation that passes `--rng-seed` to `clock` or `mugration` now fails with an unknown-argument error instead of silently ignoring the value

## Consequences

- Seeded runs are reproducible within v1, but not identical to v0, because v1 uses the Isaac64 generator and v0 uses NumPy's PCG64 ([ancestral-sample-mode-default-argmax.md](ancestral-sample-mode-default-argmax.md), [timetree-stochastic-polytomy-resolution.md](timetree-stochastic-polytomy-resolution.md))
- Tests draw random numbers only through `get_random_number_generator()`; `clippy.toml` rejects `rand::thread_rng`, `rand::rngs::StdRng`, `rand::rngs::SmallRng`, and `SeedableRng::from_entropy`
