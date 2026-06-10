# Node data JSON omits input-tree branch confidence

The timetree node data JSON (`timetree.augur-node-data.json`) and the optimize node data JSON (`optimize.augur-node-data.json`) never emit the per-node `confidence` field. Augur's `augur refine` reads `confidence` from the input tree via `getattr(n, 'confidence')` (a branch support value parsed from the input Newick/Nexus), then copies it into the node data, in both the timetree and non-timetree paths (`attributes = ['branch_length', 'confidence']`); `augur export v2` colors by `node_attrs.confidence.value`.

v1's Newick reader (`packages/treetime-io/src/nwk.rs`) does not parse internal-node support/confidence values into the node payload, so both the timetree writer (`packages/treetime/src/commands/timetree/output/augur_node_data.rs`) and the optimize writer (`packages/treetime/src/commands/optimize/augur_node_data.rs`) always set `confidence` to `None`. When the input tree carries no support values augur also omits the field, so output matches in that case. The gap appears only when the input Newick has branch support values: v1 silently drops them.

This is distinct from the auspice pseudo-bootstrap `confidence` (`1 - exp(-n_mutations)`) tracked in [Auspice JSON output missing mutations, branch confidence, and genome annotations](N-timetree-auspice-json-incomplete.md); that value is computed, whereas this one is read from the input tree.

Fixing requires the Newick/Nexus reader to parse branch support values into a node payload field, which would also benefit other commands. The `util-newick` parser already extracts support/confidence annotations; they need to be wired through `NodeFromNwk` (tracked in [kb/tickets/io-nwk-comment-dialect-parsing.md](../tickets/io-nwk-comment-dialect-parsing.md)).
