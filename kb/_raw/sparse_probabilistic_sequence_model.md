# Sparse probabilistic sequence models

In many use cases for TreeTime, branches are short and only a small minority of sites change on a given edge or in the vicinity of a node.
Calculating propagators for the majority of sites that have a negligible likelihood of changing is wasteful.
One possible approach to avoid this computational burden is to use parsimony and forgo the probabilistic models.
Another, employed by MAPLE, is to restrict probabilistic calculations to sites that are more likely to vary.
MAPLE implements this idea via Mutation Annotated Trees (MATs) and keeps track of all sites that differ from the root.
This representation ceases to be "sparse" as the tree gets deeper.
Ideally, a locally sparse representation would only keep track of sites that are variable in some radius.
This information is in principle available via the compressed representation (parsimony reconstruction) in the form of mutations on edges and ambiguous states.
Deciding what positions are to be kept variable, and what "fixed" state they are connected to outside the variable bubble, requires care.

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

For each leaf, we use the composition of the sequence at "fixed" sites from the parsimony reconstruction and instantiate a common probability distribution `[0,1,0,0]` for each state A, C, G, T from the profile map.
This `fixed` data structure tells us how likely a site has mutated away from the `fixed` state recorded in the `fixed_counts`.
For each of the sites recorded as variable in the parsimony, we record:

```
dis: [1, 0, 1, 0]
state: 'G'
```

where the `state` is the reconstructed parsimony state and the `dis` is the profile vector of the variable site.

#### Internal nodes

As we traverse the tree, the set of positions considered variable needs updating.
This is done by populating a map from positions to parsimony state from the following sources:

- Substitutions on edges to the children: each position with a substitution is added, the reference state of the substitution used as state. The state of the child is recorded as well in a way that is specific to the child.
- All positions that are variable in the children. The parsimony state of the node at these positions is the same as the child unless there is a substitution at this position along the edge to the child.

Since these variable positions are gathered from different children and edges towards them, another loop over the children collects the state of all children at these sites.
These child states are needed when the site is not variable in the child, since in this case the message used is the `fixed`-site message corresponding to the state.
For all positions that do not mutate on an edge to a child and are not variable in the child, this state has to be the parent state.

Once variable sites and the corresponding states are collected, the messages can be combined into a message that the node in question can relay to its own parent.
The `combine_messages` function calls the states in each child `reference_states` (could be renamed to `fitch_states` for clarity; it should not be called `child_states` since the same function is used in the forward pass).
The `combine_messages` function multiplies messages from different children for variable and fixed positions.
It drops variable positions when (i) all messages agree on the same parsimony state and peak probability is above a threshold. Requiring the argmax to agree with the parsimony state would be more conservative, though this condition is expected to hold in practice.

## Forward pass

The forward pass operates on the edge of the node to its parent and uses all positions that are variable in the parent, mutate on the edge, and are variable in the message to the parent as variable.
The propagated `msg_to_child` and the `msg_to_parent` are combined to yield the profile (marginal probability distribution) of variable and fixed sites at the node.

The `msg_to_child` of all children is computed by dividing out `msg_from_child` from the profile.
