# Newick output with annotations is incompatible with Newick reader

Commands that write annotated Newick trees (`optimize`, `ancestral`, `timetree`, `mugration`) produce output that the v1 Newick reader rejects. The `bio::io::newick` parser fails on any `[...]` content instead of treating it as a discardable comment per the Newick specification.

The `optimize` command writes `[&mutations="G169603A"]` annotations via `MutationCommentProvider` and `CommentProviders`, producing files like:

```
(A:0.1[&mutations="G42T"],B:0.2)root;
```

Feeding this tree back into `ancestral` (or any other command) produces a parse error at the first `[` character:

```
Error while parsing tree:  --> 1:21
  = expected SubTree
```

## Root cause

The `bio` crate (v2.3.0) Newick parser does not handle comments at all. The Newick and NEXUS specifications define `[...]` as comments that parsers should ignore, but `bio::io::newick::read` treats `[` as an unexpected token and errors.

## Affected commands

All commands writing via `write_graph_files_with_options` with non-empty `CommentProviders`:

- `optimize` (mutations) at [packages/treetime/src/commands/optimize/run.rs#L59-L66](../../packages/treetime/src/commands/optimize/run.rs#L59-L66)
- `ancestral` (mutations) at [packages/treetime/src/commands/ancestral/run.rs#L154-L176](../../packages/treetime/src/commands/ancestral/run.rs#L154-L176)
- `timetree` (mutations + dates) at [packages/treetime/src/commands/timetree/run.rs#L106-L107](../../packages/treetime/src/commands/timetree/run.rs#L106-L107)
- `mugration` (discrete states) at [packages/treetime/src/commands/mugration/run.rs#L74-L77](../../packages/treetime/src/commands/mugration/run.rs#L74-L77)

Both `.nwk` and `.nexus` output are affected since `.nexus` embeds the same annotated Newick string.

## Reproduction

```bash
./dev/docker/run ./dev/dev r treetime -- optimize --tree=data/mpox/clade-ii/1000/tree.nwk --aln=data/mpox/clade-ii/1000/aln.fasta.xz --outdir=tmp/optimize/mpox
# produces tmp/optimize/mpox/annotated_tree.nwk with [&mutations="..."] annotations

./dev/docker/run ./dev/dev r treetime -- ancestral --method-anc=marginal --tree=tmp/optimize/mpox/annotated_tree.nwk --aln=data/mpox/clade-ii/1000/aln.fasta.xz --outdir=tmp/ancestral/mpox
# fails: "expected SubTree" at position of first [&mutations=...]
```

Workaround: strip annotations with `sed 's/\[&[^]]*\]//g' input.nwk > clean.nwk`.

## Related

- [kb/issues/N-io-write-graph-files-missing-formats.md](N-io-write-graph-files-missing-formats.md)
- [kb/reports/newick-annotation-dialects.md](../reports/newick-annotation-dialects.md)
