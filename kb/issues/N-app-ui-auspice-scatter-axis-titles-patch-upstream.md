# Auspice titles scatterplot axes of branch lengths by their option labels

## Problem

Auspice 3.2.0 gives an axis of a branch length the title of that branch length: "Date", "Divergence" or "Mutations" (by the values on the axis), or the number of mutations of a gene. In the scatterplot layout it tests the labels of the axis variables (such as "Divergence" or "Mutation counts: nuc") instead of the variables (`div`, `mutationCount:nuc`), so the test never matches and the axes show the option labels. A divergence axis therefore reads "Divergence" also when its values are mutations, unlike the same branch length in the other layouts. The branch for the y axis also guesses the title from the x axis.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: the scatterplot axis titles test the axis variables (`scatterVariables.x`, `scatterVariables.y`) and guess from the scale of their own axis.

## Code

- Scatterplot axis titles: [tree/phyloTree/grid.js#L416-L424](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/tree/phyloTree/grid.js#L416-L424)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
