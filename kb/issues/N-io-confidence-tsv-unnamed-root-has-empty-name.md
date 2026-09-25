# The date interval table writes an unnamed root with an empty name

## Summary

When the input tree leaves the root unnamed, `timetree` names it `NODE_0000000` in the Auspice JSON and the other tree outputs, but the date interval table (`*.confidence.tsv`) writes the same node with an empty `name` cell. A consumer that joins the table to the tree by node name cannot find the root, the node whose date interval matters most.

## Evidence

- `fn extract_confidence_intervals()` takes the name as `names[&key].clone().unwrap_or_default()`, so a node without an input name gets the empty string [`packages/treetime/src/timetree/confidence.rs#L112`](../../packages/treetime/src/timetree/confidence.rs#L112)
- On `data/zika/86` with `--confidence --covariation`, `timetree.auspice.json` names the root `NODE_0000000`, while `timetree.confidence.tsv` has a row with an empty name and no `NODE_0000000` row; every other node name agrees between the two files

## Impact

- Joining the interval table with the Auspice tree, the node data or the Newick tree by name misses the root
- The app's result views read intervals from the Auspice JSON, so they are not affected
