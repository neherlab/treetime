# Pruning (v1-Only Command)

- [x] **Prune command** (no v0 standalone equivalent)
- [x] Tree input from `--tree`
- [x] Optional alignment input from `--aln`
- [x] Output pruned Newick and Nexus trees

## Pruning Criteria

- [x] Short branch pruning (`--prune-short` threshold)
- [x] Empty branch pruning (`--prune-empty`, requires alignment)
- [x] Name-based pruning (`--prune-nodes-list`, `--prune-nodes-list-file`)
- [x] Custom delimiters for node lists

## Validation

- [x] `--prune-empty` without `--aln` returns explicit error
- [x] `--prune-empty` loads sparse partitions and compresses sequences
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
