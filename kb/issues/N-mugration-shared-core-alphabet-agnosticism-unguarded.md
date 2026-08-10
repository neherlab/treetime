# Shared marginal/count core is alphabet-agnostic but the invariant is unguarded

## Symptom

Mugration (discrete-trait reconstruction) reuses the nucleotide reconstruction machinery: the same two-pass marginal core [packages/treetime/src/partition/marginal_core.rs](../../packages/treetime/src/partition/marginal_core.rs) and the same transition-count path (`MarginalData::count_transitions`, consumed through the `TransitionCounting` trait). Correct discrete-trait results depend on that shared core staying alphabet-agnostic. No test asserts the agnosticism contract, so a later nucleotide-specific change to the shared core (a hardcoded alphabet size, an ACGT or gap assumption, fixed 4-state indexing) could silently produce wrong discrete-trait posteriors or wrong discrete GTR estimates.

This is a coverage and hardening gap, not a demonstrated defect. The current code is agnostic; the concern is that nothing locks it that way.

## Reproduction

None. There is no failing behavior today. The item is the absence of a guard: change any shared-core computation to assume four states, and only the single 2-state discrete marginal test would catch it (and the count path would not be caught at all).

## Current state (verified, no present leak)

The shared core derives state count from the data and the model, not from a nucleotide constant:

- Marginal normalization uses profile dimensions: `dis.fill(1.0 / dis.len() as f64)` and `row.fill(1.0 / n_cols)` [packages/treetime/src/partition/marginal_core.rs](../../packages/treetime/src/partition/marginal_core.rs).
- The shared count path sizes from the model alphabet: `let n_states = self.gtr.pi.len();` with `nij` shaped `(n_states, n_states)` and `Ti` shaped `(n_states,)` [packages/treetime/src/partition/marginal_core.rs#L477-L479](../../packages/treetime/src/partition/marginal_core.rs#L477-L479).
- The root-uniform filter threshold is `1/n_states + 1e-10`, alphabet-agnostic (see [../decisions/mugration-root-state-filtering.md](../decisions/mugration-root-state-filtering.md)).
- The discrete partition supplies its own discrete `GTR`, state count `n_states()`, and leaf boundaries `one_hot_profile(index, n_states)` / `uniform_profile(n_states)` [packages/treetime/src/partition/marginal_discrete.rs](../../packages/treetime/src/partition/marginal_discrete.rs).

So the discrete path plugs into an alphabet-parameterized core through the `PartitionMarginalPasses` and `TransitionCounting` traits; no nucleotide constant is currently reachable from it.

## Impact and scope

- Severity N: no demonstrated runtime effect. The risk is a future silent-correctness regression in mugration, and in any other non-nucleotide partition, if the shared core acquires an alphabet assumption.
- The marginal recurrence has incidental regression coverage from the discrete end-to-end test; the shared count path relied on by discrete GTR estimation does not.

## Coverage gap

- `packages/treetime/src/commands/mugration/__tests__/test_discrete_marginal.rs` runs the shared marginal core on a 2-state discrete alphabet and asserts node profiles. It exercises one small state count and does not assert the agnosticism contract or a non-4 alphabet through the count path.
- The discrete transition-count delegation lacks a complete hand oracle; this narrower part is tracked in [N-mugration-count-transitions-untested.md](N-mugration-count-transitions-untested.md).

## Fix approach

- Add a regression test that drives the shared marginal and count paths on discrete alphabets with `n_states != 4` (for example 3 and 5 states) and asserts posteriors and `nij` / `Ti` against a hand oracle. This pins the invariant against an accidental nucleotide assumption and complements the count-oracle work in [N-mugration-count-transitions-untested.md](N-mugration-count-transitions-untested.md).
- Optionally add a debug assertion at the shared-core seam that profile and count dimensions equal `gtr.pi.len()`, documenting the alphabet-parameterized contract where it is relied upon.

## Related

- [N-mugration-count-transitions-untested.md](N-mugration-count-transitions-untested.md) -- narrower count-path oracle gap
- [M-discrete-missing-zero-states-inf.md](M-discrete-missing-zero-states-inf.md) -- discrete constructor admits `n_states == 0`
- [../decisions/mugration-root-state-filtering.md](../decisions/mugration-root-state-filtering.md) -- root-uniform filter, per-partition policy
- [../algo/mugration.md](../algo/mugration.md) -- shared marginal infrastructure with discrete partition
