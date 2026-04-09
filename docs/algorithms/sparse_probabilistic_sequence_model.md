# Sparse probabilistic sequence models
In many use cases for TreeTime, branches are short and only a small minority of sites change on a given edge or in the vicinity of a node.
Calculating propagators for the majority of sites that have a negligible likelihood of changing is wasteful.
One possible approach to avoid doing this computational burden is to use parsimony and forgoe the probabilistic models.
Another, employed by MAPLE, is to restrict probabilistic calculations to sites that are more likely to vary.
MAPLE implements this idea via Mutation Annotated Trees (MATs) and keep track of all sites that differ from the root.
This, however, ceases to be a "sparse" representation as the tree gets deeper.
Ideally, we'd have a locally sparse representation that only keeps track of sites that are variable in some radius.
This information is in principle available via the compressed representation (parsimony reconstruction) in form of mutations on edges and ambiguous states.
But deciding what positions are to be kept variable, and what "fixed" state they are connected to outside the variable bubble requires some care.

In the following, we assume the following is in place:
 - the sequence of the root node as constructed during the parsimony reconstruction
 - the mutations on each edge as reconstructed during the parsimony reconstruction
 - for each node the ranges of gaps and undetermined sites
 - for each tip the ambiguous sites along with the choice of the parsimony reconstruction
 - for each node the sites that were considered variable during the backward pass of the parsimony reconstruction
 - counts of the different states (e.g. A, C, G, T) in the sequence of each node

## Backward pass

### Initialization of leaves

#### Leaves
For each leave, we use the composition of the sequenced at "fixed" sites from the parsimony reconstruction and instantiate a common probability distribution $[0,1,0,0]$ for each state `A, C, G, T` from the profile map.
This `fixed` data structure tells us how likely a site has mutated away from the `fixed` state recorded in the `fixed_counts`.
For each of the sites recorded as variable in the parsimony, we record
```
 dis: [1, 0, 1, 0]
 state: 'G'
```
where the `state` is the reconstructed parsimony state and the `dis` is the profile vector of the variable site.

#### Internal nodes

As we traverse the tree, the site of positions considered variable needs updating.
This is done by populating a map from positions to parsimony state from the following sources:
 - substitutions on edges to the children: each position with a substitution is added, the reference state of the sub used as state.
   The state of the child is recorded as well in a way that is specific to the child.
 - all positions that are variable in the children. The parsimony state of the node at these positions is the same as the child unless there is a substitution at this position along the edge to the child.

Since these variable positions are gathered from different children and edges towards them, we need to run another loop over the children to collect he state of all children at these sites.
These child states are needed when the site is not variable in the child, since in this case the message used is the `fixed`-site message corresponding to the state.
For all positions that don't mutate on an edge to a child and are not variable in the child, this state has to be the parent state though.
There might be room for improvement.

Once variable sites and the corresponding states are collected, the messages can be combined into a message that the node in question can relay to its own parent.
The `combine_messages` function somewhat confusingly calls the states in each child `reference_states`; renaming this to `fitch_states` might help (it shouldn't be called `child_states` since the same function is used in the forward pass).
The `combine_messages` function multiplies messages from different children for variable and fixed positions.
Importantly, it drops variable positions when (i) all message agree on the same parsimony state and peak probability is above a threshold (NOTE: we could require the argmax to agree with the parsimony state to be more conservative -- likely always true).

## Forward pass

The forward pass operates on the edge of the node to its parent and uses all positions that are variable in the parent, mutate on the edge, and are variable in the message to the parent as variable (note sure what `parent_state.entry(*pos).or_insert_with(|| m.reff());` (the `|| x`);).
The propagated `msg_to_child` and the `msg_to_parent` are combined to yield the profile (marginal probability distribution) of variable and fixed sites at the node.

Finally, the `msg_to_child` of all children are computed by 'dividing out' the msg_from_child` out of the profile.









