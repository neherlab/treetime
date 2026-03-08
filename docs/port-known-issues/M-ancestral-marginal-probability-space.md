# Marginal reconstruction uses plain probability space

v1 dense forward pass operates in plain probability space (multiply/divide
probabilities, normalize by sum). v0 preorder operates in neg-log space (add/subtract
neg-log probabilities). Both compute mathematically equivalent operations, but the
floating-point paths differ.

The backward pass was converted to log space and matches v0.

- v1 dense forward pass: divides child message in probability space
  in [`packages/treetime/src/representation/partition/marginal_dense.rs#L213-L253`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L213-L253)
- v0 preorder: divides in log-space, multiplies back in probability space
  at [`packages/legacy/treetime/treetime/treeanc.py#L880-L917`](../../packages/legacy/treetime/treetime/treeanc.py#L880-L917)

Division by near-zero probabilities in the forward pass can amplify numerical
errors for positions where a child contributes near-zero likelihood for some
states. v0 avoids this by operating in log space (subtraction instead of
division).

Golden master tests currently compare v1 dense marginal against v0 with
tolerance 1e-6 to 1e-7, absorbing differences from the representation change.
