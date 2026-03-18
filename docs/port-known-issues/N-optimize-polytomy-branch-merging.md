# Merge branches in polytomies sharing mutations not implemented

The design document (`docs/algorithms/optimize.md:20`) describes merging branches in polytomies that share the same mutations as "a known problem in tree builders for which we have ad-hoc scripts."

## Background

Tree inference programs produce arbitrary binary resolutions of polytomies. When two branches share the same mutations (identical substitution sets), they represent redundant resolutions that should be merged. This cleaning step reduces tree builder artifacts before downstream analysis.

## Current state

No v1 implementation. No v0 formal implementation (described as "ad-hoc scripts"). The `prune` command removes short and empty branches but does not detect shared mutations between sibling branches in polytomies.
