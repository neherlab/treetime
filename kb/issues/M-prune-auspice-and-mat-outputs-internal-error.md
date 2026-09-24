# Prune fails with an internal error when writing Auspice or MAT outputs

`treetime prune` stops with an internal error when the output selection includes `auspice`, `mat-pb` or `mat-json`:

```
Error:
   0: edge_subs() called before marginal inference populated subs_ml for edge GraphEdgeKey(32)
```

The other prune outputs (`nwk`, `nexus`, `graph-json`, `dot`, `gtr`) are written without error, and `optimize` writes `auspice` and `mat-pb` from the same input. The default prune selection excludes these three formats, so only an explicit `--output-selection` hits the error.

## Reproduction

```bash
treetime prune --tree=data/flu/h3n2/20/tree.nwk --aln=data/flu/h3n2/20/aln.fasta.xz --prune-empty \
  --output-all=<dir> --output-selection=auspice
```

This is the `prune/flu/h3n2/20/auspice-mat-outputs` case of `dev/smoke`, declared an expected failure in `dev/smoke.toml`.

## Mechanism

Prune builds a Fitch partition, converts it with `into_marginal_sparse`, and seeds a `SparseReconstruction` without running marginal inference ([packages/treetime/src/prune/pipeline.rs#L54-L56](../../packages/treetime/src/prune/pipeline.rs#L54-L56)). The Auspice and MAT writers read edge substitutions through the sparse partition's `edge_subs`, which requires the maximum-likelihood substitutions (`subs_ml`) that only marginal inference fills, and errors otherwise ([packages/treetime/src/partition/marginal/sparse/partition.rs#L45](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L45)). The fix is either to provide the Fitch substitutions to these writers for prune, or to reject these formats for prune with a user-facing error.
