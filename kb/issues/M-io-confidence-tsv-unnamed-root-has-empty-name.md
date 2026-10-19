# The date interval table writes an unnamed root with an empty name

## Summary

When the input tree leaves the root unnamed, `timetree` names it `NODE_0000000` in the Auspice JSON and the other tree outputs, but the date interval table (`*.confidence.tsv`) writes the same node with an empty `name` cell. A consumer that joins the table to the tree by node name cannot find the root, whose date interval is the age of the most recent common ancestor of the samples.

## Cause

The timetree pipeline extracts the date intervals before it names the nodes:

- `fn gather_results()` calls `extract_confidence_intervals()` with `state.names` and calls `assign_node_names()` only afterwards [packages/treetime/src/timetree/pipeline.rs#L455-L479](../../packages/treetime/src/timetree/pipeline.rs#L455-L479)
- `fn extract_confidence_intervals()` takes the name as `names[&key].clone().unwrap_or_default()`, so a node without an input name gets the empty string [packages/treetime/src/timetree/confidence.rs#L100](../../packages/treetime/src/timetree/confidence.rs#L100)

Only a polytomy-resolution round that changes the topology names every node before the extraction [packages/treetime/src/timetree/round.rs#L162](../../packages/treetime/src/timetree/round.rs#L162). So the empty names occur in every run without `--resolve-polytomies` (the default), and in runs where polytomy resolution leaves the topology unchanged. The empty name applies to every node without an input name, not only the root.

## Evidence

On `data/zika/86`:

```bash
treetime timetree --tree=data/zika/86/tree.nwk --metadata=data/zika/86/metadata.tsv \
  --alignment=data/zika/86/aln.fasta.xz --time-marginal=only-final \
  --output-all=<dir> --output-confidence-tsv=<dir>/timetree.confidence.tsv
```

`timetree.auspice.json` names the root `NODE_0000000`. The interval table has one row with an empty name and no `NODE_0000000` row; every other node name agrees between the two files.

## Fix direction

Assign the node names before the interval extraction.

## Impact

- Joining the interval table with the Auspice tree, the node data or the Newick tree by name misses the root
- The app's result views read intervals from the Auspice JSON, so they are not affected
