# Graphviz node labels are not escaped

`fn print_node()` [packages/treetime-io/src/graphviz.rs#L85](../../packages/treetime-io/src/graphviz.rs#L85) writes each node as `{key} [label="({key}) {name}"]` and inserts the node name as it is. DOT strings end at an unescaped double quote, so a name that contains `"` or ends with `\` produces an invalid `.dot` file. Newick input allows such names in quoted labels, for example `('A"B':0.1,C:0.2);`.

## Expected behavior

Escape `\` and `"` in node names before writing them into a DOT string, so every `.dot` output parses whatever the node names contain.

## Validation

- A unit test writes a graph with a name that contains `"` and `\` and checks the escaped label line.
