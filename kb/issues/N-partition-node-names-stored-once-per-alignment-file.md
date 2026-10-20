# Node names are stored once per alignment file

Each alignment file of a run becomes one `struct SequenceSource` with a map from node key to `struct NodeSeqInput { name, seq }` [packages/treetime/src/seq/alignment.rs#L32-L38](../../packages/treetime/src/seq/alignment.rs#L32-L38), [packages/treetime/src/seq/alignment.rs#L122-L126](../../packages/treetime/src/seq/alignment.rs#L122-L126). The node names are the same in every file, because they come from the tree, so a run with $k$ files holds $k$ copies of every node name. The partitions read the name for error messages, for example "Leaf sequence not found" [packages/treetime/src/partition/marginal/dense/partition.rs#L53](../../packages/treetime/src/partition/marginal/dense/partition.rs#L53) and "Leaf sequence not found after alignment completion" in the Fitch pass [packages/treetime/src/partition/fitch/passes.rs#L66](../../packages/treetime/src/partition/fitch/passes.rs#L66).

The cost is one short string per node per file, small against the sequences. The second representation of the names can drift from the names of the tree.

## Fix direction

Keep node names in one place, the names of the tree, and store only sequences per source; partitions that need a name for a message take it from that one place.
