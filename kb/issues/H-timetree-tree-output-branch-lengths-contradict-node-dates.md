# Timetree tree outputs use branch lengths that contradict the node dates

The time tree written by `treetime timetree` (Newick, Nexus, MAT) has branch lengths that do not equal child date minus parent date. The node dates in the same files, in Auspice JSON, and in augur node data are the inferred dates. A reader who sums the Newick branch lengths from the root gets different dates than the ones written beside them.

## Evidence

Smoke snapshot of commit `5aeba666`, case `timetree/ebola/100/basic`:

| Leaf        | Newick length | Parent date | Leaf date | Date difference | augur `branch_length` |
| ----------- | ------------- | ----------- | --------- | --------------- | --------------------- |
| `EM_079497` | `0.0616`      | `2014.153`  | `2014.26` | `0.107`         | `0.1073`              |

The Newick lengths in this case take values close to multiples of about `0.06` years, which is one substitution divided by the clock rate. This is the pattern of a per-branch sequence-only estimate, not of a date difference.

## Mechanism

- `fn write_timetree_tree_outputs()` uses `edge.time_length` as the Newick and Nexus weight ([packages/app-output/src/timetree_tree_output.rs#L34-L35](../../packages/app-output/src/timetree_tree_output.rs#L34-L35))
- In marginal branch-length mode, `time_length` is set to `distribution.likely_time()`, the peak of the per-edge branch-length likelihood from the sequence data alone ([packages/treetime/src/timetree/inference/runner.rs#L180-L192](../../packages/treetime/src/timetree/inference/runner.rs#L180-L192)). The only other writers are the input-branch-length mode ([runner.rs#L227](../../packages/treetime/src/timetree/inference/runner.rs#L227)) and polytomy resolution
- Nothing sets `time_length` to the inferred `t_child - t_parent` after the backward and forward passes

[kb/decisions/timetree-clock-constrained-profile-propagation.md](../decisions/timetree-clock-constrained-profile-propagation.md) already records that `time_length` is the unconstrained per-edge peak. [N-io-timetree-divergence-tree-output-unimplemented.md](N-io-timetree-divergence-tree-output-unimplemented.md) assumes that the Newick and Nexus outputs contain time-based branch lengths. That assumption is correct only for the units, not for the values.

v0 writes the time tree with branch lengths derived from the inferred node dates, so the tree and its dates are consistent.

## Impact

- Wrong scientific output on the default path: every consumer of the time-tree Newick, Nexus or MAT (tree viewers, downstream dating, skyline tools) sees a tree that is not the inferred time tree
- The outputs of one run contradict each other: Newick against Auspice and augur node data

## Fix direction

Derive the tree-output branch length from the final node dates in one place, and have all tree writers read that one value. `time_length` stays an internal quantity of the inference.

## Related issues

- [M-timetree-consumers-read-unconstrained-branch-lengths.md](M-timetree-consumers-read-unconstrained-branch-lengths.md): internal consumers that read the unconstrained length
- [N-io-timetree-divergence-tree-output-unimplemented.md](N-io-timetree-divergence-tree-output-unimplemented.md)
