# Coalescent Tc re-estimation aborts on an inverted edge from round 2

`treetime timetree --coalescent-opt` fails on `data/mpox/clade-ii/1000` at refinement round 2:

```
Error:
   0: Failed to estimate the coalescent timescale
   1: Coalescent edge has child older than parent: child key=GraphNodeKey(1864),
      child=2.022958e3, parent=2.022979e3

Location: packages/treetime/src/coalescent/edge_data.rs:98
```

The child is dated about 0.021 y (roughly 7.7 days) before its parent. `--coalescent-opt
--resolve-polytomies` fails identically. Without a coalescent the same dataset completes.

## Latent, not new

The guard in
[`edge_data.rs`](../../packages/treetime/src/coalescent/edge_data.rs) has always rejected inverted
edges. It became reachable when the refinement loop started running more than one round: the
in-loop `estimate_coalescent_tc` call is gated on `i >= 2`, so before
[timetree-convergence-on-node-times.md](../decisions/timetree-convergence-on-node-times.md) it
never executed. The same inversion previously surfaced as a silent `NaN` in `log_lh_coal`.

The inversion is produced upstream of the coalescent, in the forward pass: it clamps internal
children to their parent's time but leaves observed leaf dates as given, so a leaf whose date
conflicts with the fitted clock stays where the data put it. The refinement loop's own
`commit_clock_branch_lengths` reports the same condition on this dataset — 4 branches without a
coalescent, 5 with — and commits those lengths as zero.

## Impact

`--coalescent-opt` and `--coalescent-skyline` abort on datasets containing a leaf dated before its
parent. Any tree with dating noise on the order of the clock's resolution can contain one. This
blocks the configuration that otherwise benefits most from the loop now iterating.

## Options

The two roles of the check pull in different directions, and the choice depends on
[M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md),
which is the undecided contract for what a committed node time means.

- **Tolerate in the statistic.** `collect_coalescent_edges` already skips nodes flagged as bad
  branches; extend that to inverted edges, with a count reported once rather than per edge. Cheap
  and unblocks the configuration, but silently drops those edges from the $T_c$ estimate.
- **Clamp at the source.** Project leaf times onto the topology constraint in the forward pass as
  internal nodes already are. Changes inferred dates for conflicting leaves, which is exactly the
  behavior the marginal-projection issue has not yet approved.
- **Fail, but earlier and better.** Keep the rejection and validate at the point the inversion is
  created, so the message names the leaf and its conflicting date rather than surfacing three
  layers up inside a Tc solve.

## Reproduction

```
treetime timetree --aln data/mpox/clade-ii/1000/aln.fasta.xz \
  --tree data/mpox/clade-ii/1000/tree.nwk \
  --metadata data/mpox/clade-ii/1000/metadata.tsv \
  --coalescent-opt --max-iter 6
```

## Related issues

- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md)
- [M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md](M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md)
- [M-timetree-forward-pass-skips-uncertain-leaf-refinement.md](M-timetree-forward-pass-skips-uncertain-leaf-refinement.md)
