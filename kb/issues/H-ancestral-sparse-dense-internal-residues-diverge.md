# Sparse and dense marginal reconstruction disagree on some internal-node residues

## Symptom

For the same tree and alignment, the sparse and dense marginal backends reconstruct different residues at some positions of internal nodes. Both backends place every gap and unknown character at the same positions ([kb/decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md](../decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md)), so the remaining differences are canonical states only:

| Dataset         | Internal nodes | Nodes that differ | Positions that differ |
| --------------- | -------------- | ----------------- | --------------------- |
| `data/rsv/a/20` | 18             | 6                 | 19                    |
| `data/sc2/5100` | 5099           | 467               | 658                   |

Examples on `data/rsv/a/20` (sparse, dense; positions are 0-based alignment columns): `NODE_0000002` position 5437 `C`, `G` and positions 5444, 5458, 5467, 5479, 5488 `C`, `T`; `NODE_0000016` position 6391 `A`, `G`. Examples on `data/sc2/5100`: position 29866 `A`, `T` on `NODE_0000030` to `NODE_0000034`; position 25701 `T`, `C` on `NODE_0000040`.

The differences are canonical states at positions where both backends agree on gaps and unknown characters, so the dense gap rule does not cause them.

The ambiguity report ([kb/decisions/ambiguity-report-of-plausible-states.md](../decisions/ambiguity-report-of-plausible-states.md)) shows the same divergence in the posteriors. On `data/dengue/100` the two backends differ at 92 internal-node positions without the report. With `--ambiguity-threshold=0.2`, dense writes 318 ambiguity codes and sparse 291, and at 13 more positions dense writes a code where sparse writes one state, for example `NODE_0000023` position 21 (`R` dense, `A` sparse) and position 10388 (`Y`, `T`), in each case with the same most likely state in both backends.

## Reproduction

```bash
./dev/docker/run just r treetime ancestral --method-anc=marginal --representation=sparse --tree=data/rsv/a/20/tree.nwk --alignment=data/rsv/a/20/aln.fasta.xz --output-all=tmp/residues/sparse
./dev/docker/run just r treetime ancestral --method-anc=marginal --representation=dense --tree=data/rsv/a/20/tree.nwk --alignment=data/rsv/a/20/aln.fasta.xz --output-all=tmp/residues/dense
```

Compare `ancestral.sequences.fasta` of the two runs position by position.

## Measurement on simulated data

The property test `test_prop_marginal_dense_sparse_simulated_agree_on_states` [packages/treetime/src/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.rs#L12](../../packages/treetime/src/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.rs#L12) reconstructs data from the simulation generator (`arb_simulated_case()`: coalescent trees of 6 to 20 tips, a random GTR, 300 to 1500 sites, no gaps or unknown characters) with the true GTR and the exact branch lengths. It asserts that the most likely state of every internal node is the same in both backends at every site where the two largest dense posterior probabilities differ by more than `1e-3`. It runs 64 cases and is `#[ignore]`d because of this issue.

Measurement: 1000 cases of `arb_simulated_case()` from a deterministic proptest runner compared 9950657 node sites. 3 cases (0.3%) have one differing most likely state each, so about 17% of the runs fail on their 64 new random cases. The committed regression case below also has differing states, and proptest runs it before the new cases, so every run of the ignored test fails. One more site differs among the excluded near-tie sites. Posteriors in the order `A`, `C`, `G`, `T`:

- `node_0006` position 851, dense posterior gap 0.99984: dense `[7.2e-8, 4.5e-8, 7.8e-5, 0.99992]` (`T`), sparse `[2.5e-11, 2.1e-11, 0.99999999995, 2.3e-12]` (`G`). Both backends are almost certain of different states, so this site is not a near tie. Case `CaseSpec { seed: 12108743558584770418, sampling: [SampleBatch { tau: 0.0, count: 4 }, SampleBatch { tau: 1.0, count: 4 }, SampleBatch { tau: 2.0, count: 3 }], schedule: Piecewise(PiecewiseTc { breaks: [2.5, 3.0], values: [5.0, 0.5, 1.5] }), model: RandomGtr, mu: 0.003, length: 1276 }`
- `node_root` position 1024, dense posterior gap 0.0039: dense `[0.0145, 0.0153, 0.4832, 0.4870]` (`T`), sparse `[0.0146, 0.0154, 0.4859, 0.4842]` (`G`). Case `CaseSpec { seed: 17675809652389760864, sampling: [SampleBatch { tau: 0.0, count: 4 }, SampleBatch { tau: 0.5, count: 3 }, SampleBatch { tau: 1.0, count: 4 }, SampleBatch { tau: 1.5, count: 4 }, SampleBatch { tau: 2.0, count: 3 }], schedule: Piecewise(PiecewiseTc { breaks: [2.5, 3.0], values: [5.0, 0.5, 1.5] }), model: RandomGtr, mu: 0.01, length: 1139 }`
- `node_root` position 1115, dense posterior gap 0.011: dense `[0.0011, 0.5042, 0.0016, 0.4932]` (`C`), sparse `[0.0011, 0.4975, 0.0016, 0.4999]` (`T`). Case `CaseSpec { seed: 4645121287748719626, sampling: [SampleBatch { tau: 0.0, count: 4 }, SampleBatch { tau: 1.0, count: 3 }, SampleBatch { tau: 2.0, count: 3 }], schedule: Piecewise(PiecewiseTc { breaks: [2.5, 3.0], values: [5.0, 0.5, 1.5] }), model: RandomGtr, mu: 0.003, length: 1368 }`

The two root sites move the posterior by about 0.003 to 0.007 across a near tie. The site of `node_0006` reverses a posterior of 0.9999, which a small numerical difference cannot explain. This case and the regression case below reproduce the issue without gaps or unknown characters, so they are starting points for the fix approach below.

The regression file [packages/treetime/proptest-regressions/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.txt](../../packages/treetime/proptest-regressions/ancestral/__tests__/test_prop_marginal_dense_sparse_simulated.txt) keeps a failing case of a random run. At position 928, `node_0004` and `node_0009` have a dense `G` (posterior about 0.60, `A` about 0.40) and a sparse `A` (posterior about 0.996):

```text
CaseSpec { seed: 15242536588371401067, sampling: [SampleBatch { tau: 0.0, count: 3 }, SampleBatch { tau: 0.6666666666666666, count: 3 }, SampleBatch { tau: 1.3333333333333333, count: 3 }, SampleBatch { tau: 2.0, count: 4 }], schedule: Piecewise(PiecewiseTc { breaks: [2.5, 3.0], values: [5.0, 0.5, 1.5] }), model: RandomGtr, mu: 0.003, length: 956 }
```

## Impact and scope

The sparse backend is the default. Dense and sparse are independent representations of the same model, so a differing residue means that at least one of them deviates from the marginal maximum a posteriori state. Downstream consumers see different internal sequences and mutations depending on `--representation`.

## Root cause

Not established. Leads:

- The sparse message combination demotes a variable position to a fixed one when its posterior peak exceeds `1 - EPS` ([kb/issues/M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md), lead 4). A position near this threshold can resolve to a different state than the full dense posterior
- Ties or near-ties in the posterior, where the two backends break the tie differently
- Several rsv differences cluster near positions 5437-5488 and several sc2 differences near the 3' end, where many leaves carry unknown characters

## Fix approach

Trace one differing position (for example `data/rsv/a/20`, `NODE_0000002`, 0-based position 5437) through both backends, compare the dense posterior with the sparse variable or fixed profile at that node, and identify which backend deviates from the marginal posterior. Add a sparse-versus-dense residue regression check on a small dataset once the cause is known.

The ignored test `test_dense_sparse_homoplasy_statistics_agree_with_gaps_and_unknown_characters` in [`packages/app-commands/src/commands/homoplasy/__tests__/test_dense_sparse.rs`](../../packages/app-commands/src/commands/homoplasy/__tests__/test_dense_sparse.rs) compares the `homoplasy` statistics of both backends on `data/rsv/a/20`, where the different residues move substitutions between branches: dense counts 3055 substitutions, 2145 of them on terminal branches; sparse counts 3054 and 2125. Enable the test once the backends agree.

## Related issues

The sparse backend is the default of every run ([N-representation-infer-dense-stub.md](N-representation-infer-dense-stub.md)), and the dense backend is its independent oracle. These issues show the two backends disagreeing on default settings and may share one cause. The smallest reproducing case is the tree with one unknown leaf of [H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md](H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md), so an investigation of the others starts there.

- [H-timetree-sparse-dense-node-times-differ.md](H-timetree-sparse-dense-node-times-differ.md)
- [H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md](H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md)
- [M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md)
