# Polytomy resolution merges all children indiscriminately

v0 classifies polytomy children as "stretched" (fewer mutations than clock predicts: `mutation_length < clock_length`) or "compressed" before merging. Only stretched children are merged by default (`merge_compressed=False`). Stretched branches are the ones where an intermediate ancestor could shorten the time-branch while keeping the mutation-branch short. Compressed branches (more mutations than clock predicts) are unlikely to benefit from an intermediate ancestor.

v1 merges all children without this classification, producing different resolution outcomes than v0.

v0: `packages/legacy/treetime/treetime/treetime.py:860-868`
v1: `packages/treetime/src/timetree/optimization/polytomy.rs`

## v0 reference

```python
stretched = [c for c in clade.clades if c.mutation_length < c.clock_length]
compressed = [c for c in clade.clades if c not in stretched]

if len(stretched) < 2 and merge_compressed is False:
    return 0.0

LH = merge_nodes(stretched, isall=len(stretched) == len(clade.clades))
if merge_compressed and len(compressed) > 1:
    LH += merge_nodes(compressed, isall=len(compressed) == len(clade.clades))
```

## Required changes

- Add `clock_rate` parameter to compute `clock_length = time_length * clock_rate` per child
- Classify children using `mutation_length` (from `edge.branch_length`) vs `clock_length`
- Only merge stretched children by default
- Add `merge_compressed` option
- Implement v0's `isall` stopping condition (stop at 2 when all children in group, 1 for subset)
