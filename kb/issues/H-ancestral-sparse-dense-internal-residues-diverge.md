# Sparse and dense marginal reconstruction disagree on some internal-node residues

## Symptom

For the same tree and alignment, the sparse and dense marginal backends reconstruct different residues at some positions of internal nodes. Both backends place every gap and unknown character at the same positions ([kb/decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md](../decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md)), so the remaining differences are canonical states only:

| Dataset         | Internal nodes | Nodes that differ | Positions that differ |
| --------------- | -------------- | ----------------- | --------------------- |
| `data/rsv/a/20` | 18             | 6                 | 19                    |
| `data/sc2/4500` | 5099           | 467               | 658                   |

Examples on `data/rsv/a/20` (sparse, dense): `NODE_0000002` position 5437 `C`, `G` and positions 5444, 5458, 5467, 5479, 5488 `C`, `T`; `NODE_0000016` position 6391 `A`, `G`. Examples on `data/sc2/4500`: position 29866 `A`, `T` on `NODE_0000030` to `NODE_0000034`; position 25701 `T`, `C` on `NODE_0000040`.

The differences are canonical states at positions where both backends agree on gaps and unknown characters, so the dense gap rule does not cause them.

## Reproduction

```bash
./dev/docker/run just r treetime ancestral --method-anc=marginal --dense=false --tree=data/rsv/a/20/tree.nwk --alignment=data/rsv/a/20/aln.fasta.xz --output-all=tmp/residues/sparse
./dev/docker/run just r treetime ancestral --method-anc=marginal --dense=true --tree=data/rsv/a/20/tree.nwk --alignment=data/rsv/a/20/aln.fasta.xz --output-all=tmp/residues/dense
```

Compare `ancestral.reconstructed-nuc.fasta` of the two runs position by position.

## Impact and scope

The sparse backend is the default. Dense and sparse are independent representations of the same model, so a differing residue means that at least one of them deviates from the marginal maximum a posteriori state. Downstream consumers see different internal sequences and mutations depending on `--dense`.

## Root cause

Not established. Leads:

- The sparse message combination demotes a variable position to a fixed one when its posterior peak exceeds `1 - EPS` ([kb/issues/M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md), lead 4). A position near this threshold can resolve to a different state than the full dense posterior
- Ties or near-ties in the posterior, where the two backends break the tie differently
- Several rsv differences cluster near positions 5437-5488 and several sc2 differences near the 3' end, where many leaves carry unknown characters

## Fix approach

Trace one differing position (for example `data/rsv/a/20`, `NODE_0000002`, position 5437) through both backends, compare the dense posterior with the sparse variable or fixed profile at that node, and identify which backend deviates from the marginal posterior. Add a sparse-versus-dense residue regression check on a small dataset once the cause is known.

The ignored test `test_dense_sparse_homoplasy_statistics_agree_with_gaps_and_unknown_characters` in [`packages/app-commands/src/commands/homoplasy/__tests__/test_dense_sparse.rs`](../../packages/app-commands/src/commands/homoplasy/__tests__/test_dense_sparse.rs) compares the `homoplasy` statistics of both backends on `data/rsv/a/20`, where the different residues move substitutions between branches: dense counts 3055 substitutions, 2145 of them on terminal branches; sparse counts 3054 and 2125. Enable the test once the backends agree.
