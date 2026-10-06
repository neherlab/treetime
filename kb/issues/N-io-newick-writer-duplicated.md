# Two Newick and Nexus writers

Newick and Nexus text is written by two implementations:

- `fn write_newick()` in [packages/util-newick/src/write.rs](../../packages/util-newick/src/write.rs) and `fn nexus_to_writer()` in [packages/util-newick/src/nexus.rs](../../packages/util-newick/src/nexus.rs) write a `NewickGraph`, with branch annotations, raw comments, support values and eNewick hybrid nodes. Property tests check that their output reads back unchanged
- `fn write_nwk_tree()` in [packages/treetime-io/src/nwk.rs](../../packages/treetime-io/src/nwk.rs) and `fn nex_write()` in [packages/treetime-io/src/nex.rs](../../packages/treetime-io/src/nex.rs) write the command outputs from a `TreeView`, names, branch lengths and node comments. They share only the label quoting and the comment writers of `util-newick`

Rules implemented in one writer can diverge from the other, for example the branch length format (shortest round-trip form against 3 significant digits) and the Nexus label quoting ([kb/issues/N-io-nexus-output-uses-newick-quoting.md](N-io-nexus-output-uses-newick-quoting.md)).

> [!IMPORTANT]
> **Decision required.** Make the command outputs build a `NewickGraph` and write it with `util-newick` (one owner; costs one copy of the tree per output), or move the shared rules (number format, label quoting, Nexus layout) into `util-newick` functions that both traversals call.
