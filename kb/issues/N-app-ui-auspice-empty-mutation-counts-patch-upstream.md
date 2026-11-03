# Auspice offers mutation-count branch lengths for genes without mutations

## Problem

Auspice 3.2.0 offers a mutation-count branch length for the nucleotide sequence and for each CDS of the genome map (`meta.genome_annotations`), whether or not the tree has mutations in it. A count without any mutation on the tree places every node at zero, so the tree collapses into a vertical line. TreeTime writes a CDS into the genome map for every annotated CDS of a run with amino-acid translation, also when its sequence does not change on the tree. A second tree with the genome map of the main tree but without mutations in some of its CDSs has the same problem.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: `findAvailableDistanceMeasures()` offers a mutation count only for the nucleotide sequence or a CDS with at least one valid mutation (by Auspice's `validMutation()`) on the tree.

## Code

- Mutation-count branch lengths: [util/computeDistanceMetric.ts#L110-L118](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/util/computeDistanceMetric.ts#L110-L118)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
