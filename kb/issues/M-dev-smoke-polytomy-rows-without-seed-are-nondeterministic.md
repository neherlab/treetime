# Smoke rows that resolve polytomies without a seed change between runs

`timetree --resolve-polytomies` samples the merger order at random, and only `--seed` makes the result reproducible ([kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md)). Three rows of `dev/smoke.toml` pass `--resolve-polytomies` without `--seed`:

- `timetree` `resolve-polytomies`
- `timetree` `resolve-polytomies-single-child`
- `timetree` `coal-polytomy`

Smoke compares the outputs of two snapshots byte for byte, so these rows can report `changed` without a code change. On `resolve-polytomies-single-child`, three runs of the same binary gave the root dates 1996.517, 1996.515 and 1996.848, and one smoke comparison reported 108 changed values in `timetree.augur-node-data.json`. A second comparison against the same baseline reported the case as identical.

```
treetime timetree --tree=data/smoke/flu-h3n2-20-tree-single-child-polytomy.nwk --dates=data/flu/h3n2/20/metadata.tsv \
  --aln=data/flu/h3n2/20/aln.fasta.xz --output-all=<dir> --include-leaves --resolve-polytomies
```

## Fix

Add a fixed `--seed` to each of these rows, as the sampling rows of `ancestral` and the relaxed-clock rows of `timetree` already do. The datasets of the other two rows may have no polytomy that the sampler resolves; the seed makes them reproducible either way.
