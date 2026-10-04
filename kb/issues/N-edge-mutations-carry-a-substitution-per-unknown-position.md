# Edge mutations carry a substitution per unknown position

## Symptom

The per-edge mutation lists that `ancestral` and `timetree` derive hold one substitution for every unknown position (`N`) where the other end of the edge has a residue, such as `A123N` on the edge into a sample with an `N` at position 123. Most of these entries are removed before any output is written, but they are allocated and filtered first.

On `data/sc2/2844` (4456 edges), sparse `ancestral` derives 1073122 substitutions, which `--report-ambiguous` writes to `augur-node-data`. Without that flag the same file reports 9662. The difference, about 239 substitutions per edge on average, is held in memory at the output-writing peak and then dropped.

## Mechanism

- `fn sequence_subs()` in `packages/treetime/src/seq/mutation.rs` (dense reconstructions and sampled sequences) and `fn sparse_edge_mutations()` in `packages/treetime/src/partition/marginal/sparse/mutations.rs` (sparse reconstructions) compare the two ends of an edge and skip gap positions, but not unknown positions. Every position inside a sample's unknown ranges therefore yields a substitution to `N`
- `UnknownMutationFilter` in `packages/app-output/src/mutation_filter.rs` removes these entries when it writes outputs, unless `--report-ambiguous` is set, in which case they are written as they are
- `fn bridge_edge()` in the same file uses them: a substitution to `N` records the hidden residue, and a later substitution from `N` on a descendant edge reports the change from that hidden residue, so a change across a masked stretch appears on the edge where the residue is observed again

## Reproduction

Run `treetime ancestral --method-anc=marginal --dense=false` on `data/sc2/2844`, once with `--report-ambiguous` and once without, and compare the number of substitutions in `augur-node-data`:

```bash
jq '[.nodes[] | (.muts // [] | length)] | add' ancestral.augur-node-data.json
```

The difference is the number of substitutions to or from `N` that the default output drops.

## Fix direction

Represent unknown stretches as ranges instead of one substitution per position, or skip them during derivation and give `fn bridge_edge()` the hidden residues from the node states directly. Either way, the bridged changes and the `--report-ambiguous` output must stay byte-identical. Verify with the output-consistency property test in `packages/treetime/src/ancestral/__tests__/test_output_consistency_prop.rs`, the sparse mutation property test in `packages/treetime/src/ancestral/__tests__/test_sparse_mutations_prop.rs`, and `just smoke`.
