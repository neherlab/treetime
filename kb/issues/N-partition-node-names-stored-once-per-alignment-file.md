# Node names are stored once per alignment file

Each alignment file of a run becomes one `struct SequenceSource` with a map from node key to `struct NodeSeqInput { name, seq }` [packages/treetime/src/seq/alignment.rs#L28-L34](../../packages/treetime/src/seq/alignment.rs#L28-L34), [packages/treetime/src/seq/alignment.rs#L110-L113](../../packages/treetime/src/seq/alignment.rs#L110-L113). The node names are the same in every file, because they are copies of the names of the tree (`struct NodeNames`), so a run with $k$ files holds $k$ copies of every node name. The partitions read the name for error messages, for example "Leaf sequence not found" [packages/treetime/src/partition/marginal/dense/partition.rs#L50](../../packages/treetime/src/partition/marginal/dense/partition.rs#L50) and "Leaf sequence not found after alignment completion" in the Fitch pass [packages/treetime/src/partition/fitch/passes.rs#L69](../../packages/treetime/src/partition/fitch/passes.rs#L69).

The cost is one short string per node per file, small against the sequences. The second representation of the names can drift from the names of the tree.

## Fix direction

Keep node names in one place, the names of the tree, and store only sequences per source; partitions that need a name for a message take it from that one place.
