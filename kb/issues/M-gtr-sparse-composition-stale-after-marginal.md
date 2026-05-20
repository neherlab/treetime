# Sparse GTR inference mixes MAP mutations with Fitch-era compositions

After `PartitionBranchOps` promotion, `get_mutation_counts_sparse()` derives `nij` from MAP-state mutations via `edge_subs_from_graph()`, but `Ti` base counts and `root_state` still read `SparseNodePartition.seq.composition` which is populated by Fitch compression and not updated by marginal inference.

## Impact

The `Ti` (time-in-state) base accumulates `branch_length * composition[nuc]` per edge, then adjusts for each mutation. The base counts character frequencies from Fitch compression. After marginal inference, variable-site assignments can change, making the composition stale. The adjustment terms (from MAP mutations) are correct, but the base terms are not.

Practical magnitude is small: only variable sites differ between Fitch and marginal, and these are a small fraction of total positions. The `root_state` prior is similarly affected but serves as a soft prior on equilibrium frequencies where approximate values are acceptable.

## Root cause

`SparseSeqInfo.composition` is set during Fitch compression ([packages/treetime/src/ancestral/fitch.rs](../../packages/treetime/src/ancestral/fitch.rs)) and never updated. `reconstruct_node_sequence()` rewrites `seq.sequence` but not `seq.composition` ([packages/treetime/src/partition/marginal_sparse.rs](../../packages/treetime/src/partition/marginal_sparse.rs)).

## Fix

Recompute `seq.composition` from `seq.sequence` after marginal reconstruction, or derive `Ti`/`root_state` from the same MAP state that `edge_subs_from_graph()` uses.

## v0 comparison

v0 overwrites `node.cseq` (compressed sequence) after each marginal pass. The `mutations` property dynamically compares current `node.up.cseq` vs `node.cseq`, so composition data always reflects the current reconstruction state. The stale-composition problem does not exist in v0.

## Test co-dependency

Fixing this issue requires updating tests that encode the stale behavior:

- `gtr/infer_gtr/__tests__/test_contract.rs`: `test_root_state_sparse()` asserts Fitch-consensus counts after `update_marginal`
- `gtr/infer_gtr/__tests__/test_fitch.rs`: `test_get_mutation_counts_sparse()` hard-codes `root_state` without deriving from reconstructed MAP root sequence

Both tests will fail when composition is refreshed and must be rewritten to assert post-marginal state.

## Related

- [M-gtr-per-site-rate-variation](M-gtr-per-site-rate-variation.md) - per-site rate variation, another GTR inference gap
- [docs/reports/optimization-methods/](../reports/optimization-methods/) - GTR inference comparison across tools

## Scientific background

GTR parameter estimation from phylogenetic data uses the EM framework <a id="cite-1"></a>[Holmes and Rubin 2002](https://doi.org/10.1006/jmbi.2002.5405) [[1](#ref-1)]. The E-step computes expected sufficient statistics (dwell times $T_i$ and transition counts $n_{ij}$) conditioned on observed data and current parameters. The M-step re-estimates rate matrix parameters from these statistics in closed form.

TreeTime uses an ECM variant with MAP hard assignment: ancestral states are set to their most probable value rather than integrated over the full posterior. This is a generalized EM that preserves monotone likelihood increase but introduces bias proportional to ancestral state uncertainty <a id="cite-2"></a>[Meng and Rubin 1993](https://doi.org/10.1093/biomet/80.2.267) [[2](#ref-2)].

The stale composition issue means the E-step uses mixed provenance: $n_{ij}$ from current MAP states (correct) but $T_i$ from Fitch-era compositions (stale). <a id="cite-3"></a>[Puller, Sagulenko, and Neher 2020](https://doi.org/10.1093/ve/veaa066) [[3](#ref-3)] show that ignoring site-specific preferences in TreeTime's GTR inference underestimates branch lengths, with magnitude comparable to ignoring rate variation. The stale composition is a milder version of the same problem: base counts from a previous reconstruction stage fed into the current GTR inference iteration.

For typical viral datasets with low diversity, variable sites are a small fraction of total positions, bounding the discrepancy between Fitch and marginal compositions. The bias grows with sequence divergence and the proportion of ambiguous sites.

## References

1. <a id="ref-1"></a> Holmes, Ian, and Gregory M. Rubin. 2002. "An Expectation Maximization Algorithm for Training Hidden Substitution Models." _Journal of Molecular Biology_ 317(5):753-764. https://doi.org/10.1006/jmbi.2002.5405 [↩](#cite-1)
2. <a id="ref-2"></a> Meng, Xiao-Li, and Donald B. Rubin. 1993. "Maximum Likelihood Estimation via the ECM Algorithm: A General Framework." _Biometrika_ 80(2):267-278. https://doi.org/10.1093/biomet/80.2.267 [↩](#cite-2)
3. <a id="ref-3"></a> Puller, Vadim, Pavel Sagulenko, and Richard A. Neher. 2020. "Efficient Inference, Potential, and Limitations of Site-Specific Substitution Models." _Virus Evolution_ 6(2):veaa066. https://doi.org/10.1093/ve/veaa066 [↩](#cite-3)
