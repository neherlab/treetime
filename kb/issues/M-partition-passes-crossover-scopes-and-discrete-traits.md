# Partition passes do not support crossover scopes or discrete traits on segment graphs

Fitch parsimony and the dense and sparse marginal passes read edge scopes ([kb/algo/ancestral.md](../algo/ancestral.md)) for edges that carry a segment entirely or not at all. Two cases remain:

- **Crossover scopes**: an edge that carries part of a segment (site ranges) is rejected by `Inheritance::require_one_parent()` ([packages/treetime-graph/src/scope.rs](../../packages/treetime-graph/src/scope.rs)). Every pass then needs, per site, the parent that carries it: the backward step can treat sites outside a child's scope as sites without data below that child, and the forward step must assemble the message from above site by site
- **Discrete traits**: mugration builds its partition with `Inheritance::ALL`, so a trait cannot follow one segment of a graph of segment trees; a trait needs a rule for the segment it follows

## Validation

- Crossover: a toy graph with one recombinant node whose two parents split a segment, against two runs on the trees of the two intervals
- Tree inputs: outputs unchanged
