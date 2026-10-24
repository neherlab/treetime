# Date-free rooting sets the split of the root edge without dates in dated runs

In `clock` and `timetree`, a date-free method (`min-dev`, `--outgroup`, `midpoint`) chooses the root from the branch lengths alone, and the root stays for the rest of the run ([kb/decisions/reroot-date-free-methods-root-once.md](../decisions/reroot-date-free-methods-root-once.md)). The method sets two things: the edge that holds the root, and the split of that edge into the two branches below the new root. The dates decide neither.

The edge is the choice the user asked for. The split is a by-product of the method: `midpoint` puts the root at the middle of the longest tip-to-tip path, `--outgroup` at the middle of the branch above the named node, and `min-dev` at the minimum of the variance of the root-to-tip distances. With `--branch-length-mode=input`, the time inference turns each branch length into a likelihood of the branch's time length (`fn input_branch_model()`, `packages/treetime/src/timetree/inference/input_branch_length.rs`), so the split goes into the likelihoods of the root's two child edges, and it can move the root date.

`least-squares` places the root on the same kind of edge with a split that fits the dates: the root-to-tip regression puts the root where the dated tips lie closest to a line.

> [!IMPORTANT]
> **Investigation required.** It is not measured whether the date-free split biases the root date. Evidence to collect: on `data/lassa/L/20`, run `timetree --reroot=midpoint --branch-length-mode=input` once as it is, and once with the root moved by hand to the least-squares split of the same edge (the split that `least-squares` gives when it chooses that edge, or the least-squares optimum on that edge), and compare the root date and its confidence interval of the two runs.

## Locations

- `fn reroot_date_free()` in `packages/treetime/src/reroot/date_free.rs`: the date-free root search and topology change
- `fn midpoint_root()` in `packages/treetime/src/reroot/midpoint.rs` and `fn named_node_point()` in `packages/treetime/src/reroot/named_node.rs`: the splits of `midpoint` and `--outgroup`
