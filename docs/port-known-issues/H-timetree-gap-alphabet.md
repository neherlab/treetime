# Gap character not handled in alphabet

**Status: resolved.**

The alphabet profile map now always includes the gap character `-` with a uniform profile (same as unknown `N`), matching v0's `nuc_nogap` behavior. The `treat_gap_as_unknown` parameter was removed as it was the sole consumer of this flag and gap-as-unknown is the only valid semantic for a 4-state (A/C/G/T) nucleotide alphabet.

The ebola_20 dataset no longer crashes during marginal dense sequence processing. The [marginal dense runner test](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L32) remains commented out due to a separate golden master node key mismatch (v0 captures 11 internal nodes, v1 rerooting produces 19).
