# Optimization

The optimize command adjusts branch lengths of an input tree to maximize the likelihood. Optimal branch lengths depend on the input alignments and the substitution models.

## Inputs

- Tree
- Alignment(s) (corresponding to partitions)
- Substitution models (standard or inferred)

For more complex inputs (multiple alignments and models along with discrete characters) a config file format may be needed.

## Algorithm

The likelihood is computed for each partition. Once all messages are in place, the optimal length of each branch is calculated using either a grid search or Newton's method from the message "to parent" and "to child" on each branch. Once all branch lengths are updated, all messages are recomputed and the procedure iterated until (approximate) convergence. The output is a new tree.

## Initial conditions

The starting point for the Newton iteration must be in a region that leads to convergence. The substitution count divided by alignment length is a good estimate in most cases. The input tree will usually provide a fine starting point already, so this step could be made optional. Only when the input tree does not have approximate branch lengths is the initial guess needed.

## Post-optimization cleanup

- Pruning branches of zero length
- Scanning for shared substitutions in children of a polytomy and introducing a new internal node with the individuals that share a substitution as children. Initial branch length for this new internal node would be the shared substitution count divided by alignment length.

## Indel contribution

An open question is whether to extend optimization to indels. Most phylogenetic software ignores indels, but the unified optimization framework allows adding an indel contribution to the likelihood. A branch with no substitutions but an indel should have a finite length, requiring indels to contribute to the likelihood.
