# Sparse ancestral outputs write the Fitch root with marginal mutations

In sparse marginal reconstruction (the default of `treetime ancestral`), the root sequence written to the tree outputs is the Fitch parsimony root, but the mutations on the branches come from the marginal reconstruction. When the two roots differ, the written root plus the written mutations does not reproduce the reconstructed sequences.

## Reproduction

`t1.nwk`:

```
((A1:0.01,A2:0.01)X:1.0,(C1:0.01,C2:0.01)Y:0.001)root;
```

`t1.fasta`: `A1` and `A2` are `AAAAGGGGAA`, `C1` and `C2` are `CAAAGGGGAT`.

```bash
treetime ancestral --tree t1.nwk --alignment t1.fasta --model jc69 --output-all out --output-selection all
```

| Output                                       | Root sequence |
| -------------------------------------------- | ------------- |
| `ancestral.auspice.json` `root_sequence`     | `AAAAGGGGAA`  |
| `ancestral.augur-node-data.json` `reference` | `CAAAGGGGAT`  |

Node `X` carries mutations `C1A` and `T10A`, and node `Y` carries none. These mutations fit the augur reference. Applied to the Auspice root, they start from the wrong states: `Y` and its tips get `AAAAGGGGAA` instead of `CAAAGGGGAT`.

With `--dense true`, both outputs write `CAAAGGGGAT`.

## Mechanism

- `PartitionFitch::into_marginal_sparse()` copies the Fitch root into `PartitionMarginalSparse.root_sequence` ([packages/treetime/src/partition/fitch/partition.rs#L31](../../packages/treetime/src/partition/fitch/partition.rs#L31)). The marginal passes do not update this field
- `MarginalReconstruction::root_sequence()` returns that field for a sparse reconstruction ([packages/treetime/src/partition/marginal/reconstruction.rs#L148](../../packages/treetime/src/partition/marginal/reconstruction.rs#L148), [packages/treetime/src/partition/marginal/sparse/partition.rs#L53-L55](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L53-L55)). The tree writers (Auspice, MAT reference) use `root_sequence()`
- `augur_root_sequence()` and `edge_subs()` read the marginal node states, which is why augur node data is consistent

## Impact

- Wrong reconstructed sequences for every consumer that rebuilds node sequences from the root and the branch mutations (Auspice, MAT-based tools)
- Sparse is the default mode, and the Fitch and marginal roots differ whenever parsimony and likelihood disagree at the root, for example with a long branch on one side of the root

## Fix direction

Return the marginal root sequence from the sparse reconstruction, so that the root and the mutations come from one source. A test that rebuilds every node sequence from `root_sequence()` plus `edge_mutations()` and compares with `augur_node_sequence()` covers the dense, sparse and Fitch paths.

## Related issues

- [M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md)
- [M-ancestral-sparse-root-invariance.md](M-ancestral-sparse-root-invariance.md)
