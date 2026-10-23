# Sparse reroot drops the unknown ranges of the new root

After a reroot of a sparse reconstruction, the new root sequence carries the states of the old root at sites where the new root is unknown (`N`). A node whose whole clade is unknown at a site changes its reconstructed state from `N` to the state of the old root when it becomes the root, so the reroot changes the sequence of a node that it only moves.

## Reproduction

Tree `((A:0.1,B:0.1)X:0.1,C:0.1,D:0.1)root;`, sequences `A=NCG`, `B=NCG`, `C=ACG`, `D=ACT`, sparse Fitch reconstruction, reroot through `apply_root_at_edge()` on the edge root->X at split 1.0 (the new root is X, the old root keeps three children and stays):

| Node | Before the reroot | After the reroot |
| ---- | ----------------- | ---------------- |
| X (parsimony sequence) | `NCG` | `ACG` |
| X (marginal `node_sequence()`) | `NCG` | `ACG` |
| root | `ACG` | `ACG` |

No edge carries an event at the first site: Fitch gives the old root `A` and records the unknown range of X in its node observation, so the root->X edge is empty there.

A reroot moves no event between nodes, so the sequence of X must stay `NCG`. The representation carries `N` in a root sequence: with every leaf `N` at the first site, the Fitch root sequence is `NCG` and the root observation has the unknown range `(0, 1)`.

## Cause

- `fn derive_root_sequence()` [packages/treetime/src/partition/marginal/sparse/reroot.rs#L145](../../packages/treetime/src/partition/marginal/sparse/reroot.rs#L145) applies only the substitutions and indels of the inverted path edges (`fn apply_edge_to_sequence()`, same file) to the old root sequence. It never applies the unknown ranges of the new root's observation
- `fn parsimony_seq()` [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L12-L37](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L12-L37), which builds every other node sequence in the forward pass, fills the node's unknown ranges with `N` after the edge events
- The sparse forward pass seeds the root state from `root_sequence` [packages/treetime/src/partition/marginal/sparse/forward.rs#L97-L101](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L97-L101), so the marginal output of the new root inherits the missing `N`
- A node split from an edge for the reroot takes its observation from the same root sequence [packages/treetime/src/partition/marginal/sparse/reroot.rs#L40-L46](../../packages/treetime/src/partition/marginal/sparse/reroot.rs#L40-L46)

## Fix direction

Fill the unknown ranges of the new root's observation into the derived root sequence after the edge events, as `parsimony_seq()` does for every other node. The invariant test `test_reroot_sparse_reconstruction_keeps_the_sequence_of_every_surviving_node` (`packages/treetime/src/reroot/__tests__/test_reroot_sparse_reconstruction.rs`) checks the parsimony sequences of every surviving node across all reroot targets; a case with a clade that is unknown at a site covers this defect.

> [!IMPORTANT]
> **Investigation required.** The reproduction ran once as a temporary test and is not committed. Before the fix, confirm the effect on a CLI run (`treetime optimize` or `timetree` with `--representation sparse` and a reroot) with a clade that is unknown at some sites, and check whether the marginal likelihood changes as well as the node sequence.
