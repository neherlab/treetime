# Optimization

the optimize command serves to adjust branch length of an input tree to maximize the likelihood. These optimal branch length depend on the input alignments and the substitution models. inputs are thus

- Tree
- alignment(s) (corresponding to partitions)
- substitutions models (standard or inferred)

For more complex inputs (multiple alignments and models along with discrete characters we might need a config file format).

The likelihood is computed for each partition. once all messages are in place the optimal length of each branch is calculated using either a grid search or Newton's method from the message "to parent" and "to child" on each branch. Once all branch lengths are update, all messages are recomputed and the procedure iterated until (approximate) convergence. The output is a new tree.

The algorithm for this is mostly in place. But some tweaks seem needed to detect when zero length is the optimal branch length (the current function doesn't work well). the run function of this command is a prototype at the moment.

From my observations optimize has a tendency to oscillate slightly which impedes convergence, the initialization is off when there are gaps or unknown characters, and the iteration is slow for dense (which comes down to ancestral being compute heavy).

Potential additional useful things:

prune zero length branches
merge branches in polytomies that share mutations (a known problem in tree builders for which we have ad-hoc scripts)

The previous paragraph on this didn't talk about initial conditions. They really only serve the purpose to start the Newton iteration from a place that will likely lead to convergence. #subs/length is a good estimate in most cases. Furthermore, the input tree will usually be a fine starting point already, so we could make this step optional. Only when the input tree doesn't have approximate branch length we need the initial guess.

More added value of the optimize command would be

- pruning branches of zero length.
- scanning for shared subs in children of a polytomy and introducing a new internal node with the individuals that share a sub as children. Initial branch length for this new internal node would be #shared subs/length.

An outstanding question is whether we can extend this to indels. we, as most other phylo software, ignore indels. But the unified optimization should allow to add an indel contribution to the LH. A branch with no subs but an indel should have a finite length. But for this to happen we need to add it somehow to the LH.
