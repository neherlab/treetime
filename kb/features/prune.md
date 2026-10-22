# Pruning (v1-Only Command)

- [x] **Prune command** (no v0 standalone equivalent)
- [x] Tree input from `--tree`
- [x] Optional alignment input from `--alignment`, one or several files
- [x] Tips without a sequence in a file are filled as unknown; the run stops when a file misses more than a third of the tips
- [x] Output pruned Newick and Nexus trees

## Pruning Criteria

- [x] Short branch pruning (`--collapse-short` threshold)
- [x] Empty branch pruning (`--collapse-empty`, requires alignment; with several files a branch is empty only when no file has a substitution on it)
- [x] Name-based pruning (`--remove`, `--remove-file`)

## Validation

- [x] `--collapse-empty` without `--alignment` returns explicit error
- [x] A given alignment loads sparse partitions and compresses sequences; the tree outputs carry its mutations
- [x] `--gap-fill` replaces terminal gaps with the ambiguous character by default, as in the other sequence commands
- [x] Mutation counts come from sparse edge substitutions

## Pruning Workflow

- [x] Two-pass strategy
  - [x] Internal-edge collapse pass
  - [x] Graph rebuild
  - [x] Leaf-removal pass
  - [x] Graph rebuild
- [x] Internal edge collapsed when: shorter than threshold, empty of substitutions, or target node name selected
- [x] Leaf pass removes only selected leaf names
- [x] Recursive parent collapse (childless parent removal)
- [x] Edge collapse with data merging (branch length sum, mutation union)
