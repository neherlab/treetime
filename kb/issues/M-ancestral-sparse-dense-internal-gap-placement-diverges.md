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

Not established. The two backends resolve indels through different code paths (`fn resolve_indels_forward()` / `resolve_indels_backward()` for dense; the sparse indel track carried on edges), and they diverge on which internal positions are deleted. A focused comparison against v0 is needed to determine which backend, if either, matches the reference and where the divergence originates.

## Fix approach

Trace one diverging internal position through both backends and against v0 (`packages/legacy/treetime/treetime/treeanc.py`), identify the correct deletion set, and reconcile the sparse indel reconstruction with it. Add a sparse-vs-dense internal-gap regression check on a small dataset once the correct behavior is known.
