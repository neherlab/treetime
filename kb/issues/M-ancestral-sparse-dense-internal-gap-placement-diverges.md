# Sparse and dense marginal reconstruction place internal-node gaps differently

## Symptom

For the same tree and alignment, the sparse and dense marginal ancestral backends disagree on where gaps (deletions) fall in internal-node sequences. On `data/rsv/a/20` with `--method-anc=marginal`, comparing the two backends' node-data JSON finds 47 internal-node positions where one backend reports a gap and the other reports a residue. The disagreement is on internal edges and is present without any tip-reconstruction flag, so it is independent of the tip reconstruction and imputation behavior.

## Reproduction

```
./dev/docker/run ./dev/dev r treetime -- ancestral --method-anc=marginal --tree=data/rsv/a/20/tree.nwk --aln=data/rsv/a/20/aln.fasta.xz --output-all=tmp/gap/sparse
./dev/docker/run ./dev/dev r treetime -- ancestral --method-anc=marginal --dense=true --tree=data/rsv/a/20/tree.nwk --aln=data/rsv/a/20/aln.fasta.xz --output-all=tmp/gap/dense
```

Then, for each shared node, count positions where exactly one backend has `-`. The internal nodes accumulate 47 such positions on this dataset.

## Impact and scope

The sparse backend is the default. Its internal-node sequences differ from the dense backend's (taken as the reference) in gap placement, so a downstream consumer of the reconstructed internal sequences sees deletions at different sites depending on the backend. Scope is indel/deletion reconstruction on internal nodes; canonical substitution states and the marginal log-likelihood are not implicated by this comparison.

## Root cause

Partly established. One contributor is fixed: sparse reconstruction re-applied only the deletions recorded on a node's *own* parent edge, so a node whose gap was inherited from further up had a residue written back by the posterior. See [The sparse stored sequence is the parsimony chain](../decisions/ancestral-marginal-sparse-parsimony-chain.md). On `data/sc2/4500` position 28369 that accounted for ~4470 internal nodes reported as a residue where dense reports a gap; sparse now matches dense there.

A residual divergence remains and is what this issue now tracks. On the same dataset, sparse reports the unknown character `N` at 25 internal nodes (position 28369) and 72 (position 23009) where dense reports a gap or a residue. That residual predates the fix above and is unchanged by it, so it has a separate cause - most likely in how the two backends classify a position as `non_char` versus deleted, rather than in the reconstruction step. The two backends resolve indels through different code paths (`fn resolve_indels_forward()` / `resolve_indels_backward()` for dense; the sparse indel track carried on edges). A focused comparison against v0 is still needed to determine which backend matches the reference.

## Fix approach

Trace one diverging internal position through both backends and against v0 (`packages/legacy/treetime/treetime/treeanc.py`), identify the correct deletion set, and reconcile the sparse indel reconstruction with it. Add a sparse-vs-dense internal-gap regression check on a small dataset once the correct behavior is known.
