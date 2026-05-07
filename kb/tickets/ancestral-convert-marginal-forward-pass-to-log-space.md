# Convert marginal forward pass to log-space arithmetic

v1 dense forward pass operates in plain probability space (multiply/divide
probabilities, normalize by sum). v0 preorder operates in neg-log space (add/subtract
neg-log probabilities). Both compute mathematically equivalent operations, but the
floating-point paths differ.

The backward pass uses log-space arithmetic in both dense (`normalize_from_log()`) and sparse (`softmax_with_log_norm()` in `combine_messages()`) modes.

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

## Related issues

- Source: [M-ancestral-marginal-probability-space.md](../issues/M-ancestral-marginal-probability-space.md) -- delete after full resolution
