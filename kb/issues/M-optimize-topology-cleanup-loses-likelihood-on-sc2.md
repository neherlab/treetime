# Topology cleanup of `optimize` returns a less likely tree than `rust` on sc2/2844

On `sc2/2844` with default settings, the final tree of `treetime optimize` re-optimized in dense mode with a fixed model (`--model=jc69 --representation=dense --no-collapse-short-branches`) has a log-likelihood of -143227, while the final tree of the previous version (`rust` before the zero-length collapse of every internal edge) reaches -143105 under the same re-optimization: 122 units less. On `flu/h3n2/500` the new tree is 0.05 units better and has the 205 clades of TreeTime v0.

Clade comparison of the two `sc2/2844` trees (inferred model): 823 internal nodes in the previous tree, 780 in the new one, 763 shared. Of the 60 clades only in the previous tree, 28 have a clade within 1 to 2 tips in the new tree (tips moved by merge or hoist), 30 have a dense stem above 0.5 mutations (up to 2.26, 35.4 in total), and 41 come from the input tree.

The previous version kept every zero-length edge that carries a Fitch substitution or an indel. The new version collapses an internal edge as soon as its optimized length is 0, also in the first iteration, when the other branch lengths are still the input values. Deferring the collapse of edges with Fitch substitutions to the collapse after the loop does not help: that tree is 393 units below the previous one.

> [!IMPORTANT]
> **Investigation required.** Find the edges that cost the likelihood: put each of the 60 clades back into the new tree, re-optimize in dense mode, and record the gain. For the edges with a dense stem above 0.5 mutations, compare their sparse and dense optimized lengths on the same tree, to tell an early collapse from a sparse objective that sets these edges to 0. The dense re-optimization of the input tree itself stops below the previous tree (inferred model), so check that 10 iterations reach the optimum before comparing trees.
