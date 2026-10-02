# Input branch-length mode drops the time evidence of edges without an input length

## Problem

With `--branch-length-mode=input`, the branch-length likelihood of an edge comes only from its input branch length. `fn create_branch_distributions_input_mode` ([packages/treetime/src/timetree/inference/runner.rs](../../packages/treetime/src/timetree/inference/runner.rs)) gives an edge whose input length is absent neither a distribution nor a time length. Two consequences follow in the time inference:

- **No fallback length**: an edge without an input branch length carries no evidence at all. Marginal mode uses one mutation as the fallback length for such an edge; input mode has no fallback. The length is absent when the input Newick omits it for a branch
- **Non-bad node without a message**: `fn send_backward_message` ([packages/treetime/src/timetree/inference/backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs)) sends no message when the edge has no branch-length distribution or the node has no subtree distribution. An internal node that is not a bad branch, because one of its children carries a date, then sends nothing to its parent when its own parent edge or every child edge lacks an input length. The parent treats the subtree as if it carried no date, although the bad-branch flags say it does

## Reachability

Both cases are latent. Input mode aborts on ordinary trees before any refinement round with `Cannot divide point by point` ([H-timetree-input-branch-lengths-abort-on-point-division.md](H-timetree-input-branch-lengths-abort-on-point-division.md)), so no input-mode run reaches them today. They become reachable when that issue is fixed.

## Open question

Decide the input-mode evidence of an edge without an input length: a fallback length (for example one mutation, as marginal mode), an uninformative branch likelihood that still forwards the subtree distribution, or an error that asks for a complete set of branch lengths. The choice also decides whether a non-bad node always sends a message to its parent.
