# Timetree polytomy resolution fails when the tree has a single-child node

`treetime timetree --resolve-polytomies` aborts with an internal error when the tree has a polytomy and an internal node with one child. The same tree without the single-child node runs to completion. The date of the single-child node does not matter: the error appears with and without a date for it.

```
0: When running round 1
1: Partition contains stale topology entries while indexing a pass. This is an internal error. Please report it to developers.
Location: packages/treetime-graph/src/pass.rs:222
```

## Reproduction

```bash
d=data/flu/h3n2/20
mkdir -p tmp/single-child

# Single-child node U above the clade of A/Oregon/15/2009 and A/Hong_Kong/H090_695_V10/2009
sed -E 's#\((A/Oregon/15/2009[^,]*,A/Hong_Kong/H090_695_V10/2009[^)]*)\)0\.977:0\.00426#((\1)0.977:0.002)U:0.00226#' \
  $d/tree.nwk > tmp/single-child/tree_u.nwk

# Polytomy: collapse the internal node with 15 tips
./dev/docker/python python -c '
from Bio import Phylo
t = Phylo.read("tmp/single-child/tree_u.nwk", "newick")
t.collapse(next(n for n in t.get_nonterminals() if n is not t.root and len(n.clades) == 2 and n.count_terminals() == 15))
Phylo.write(t, "tmp/single-child/tree_u_poly.nwk", "newick", format_branch_length="%.6f")'

# Fails with the internal error
./dev/docker/run just r treetime timetree --resolve-polytomies --tree=tmp/single-child/tree_u_poly.nwk \
  --metadata=$d/metadata.tsv --alignment=$d/aln.fasta.xz --output-all=tmp/single-child/u-poly

# Runs to completion: the same tree with U collapsed
./dev/docker/python python -c '
from Bio import Phylo
t = Phylo.read("tmp/single-child/tree_u_poly.nwk", "newick")
t.collapse(next(n for n in t.find_clades() if n.name == "U"))
Phylo.write(t, "tmp/single-child/tree_poly.nwk", "newick", format_branch_length="%.6f")'
./dev/docker/run just r treetime timetree --resolve-polytomies --tree=tmp/single-child/tree_poly.nwk \
  --metadata=$d/metadata.tsv --alignment=$d/aln.fasta.xz --output-all=tmp/single-child/poly
```

Without a polytomy (`tree_u.nwk`), polytomy resolution returns early ([packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L26-L30](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L26-L30)) and the run completes with `U` in the output.

## Analysis

After the polytomies are resolved, `remove_single_child_nodes()` removes every node with one parent and one child from the graph ([packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L212-L236](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L212-L236)), including single-child nodes of the input tree. The next pass fails the check that every node and edge key of a partition is present in the graph ([packages/treetime-graph/src/pass.rs#L219-L223](../../packages/treetime-graph/src/pass.rs#L219-L223)). This points to partition entries that remain for the removed node and its edges. The reconciliation of partitions after the topology change was not traced.

## Impact

Any input tree with a single-child internal node, such as a subtree cut from a larger tree or a tree with a sampled ancestor, cannot be used with `--resolve-polytomies` when it has a polytomy after rerooting.

## Related issues

- [M-timetree-polytomy-cleanup-drops-dated-single-child-node.md](M-timetree-polytomy-cleanup-drops-dated-single-child-node.md)
