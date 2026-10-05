# Marginal reconstruction uses plain probability space

v1 dense forward pass operates in plain probability space (multiply/divide probabilities, normalize by sum). v0 preorder operates in neg-log space (add/subtract neg-log probabilities). Both compute mathematically equivalent operations, but the floating-point paths differ.

The backward pass uses log-space arithmetic in both dense (`normalize_from_log()`) and sparse (`softmax_with_log_norm()` in `combine_messages()`) modes.

- v1 dense forward pass: `fn indexed_node_forward()` divides the parent profile by the child message in probability space
  at [`packages/treetime/src/partition/marginal/shared/pass.rs#L247-L248`](../../packages/treetime/src/partition/marginal/shared/pass.rs#L247-L248), through `fn divide_out()`, which floors the divisor at `f64::MIN_POSITIVE`
  in [`packages/treetime/src/partition/marginal/shared/normalize.rs#L56-L60`](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L56-L60)
- v0 preorder: divides in log-space, multiplies back in probability space
  at [`packages/legacy/treetime/treetime/treeanc.py#L880-L917`](../../packages/legacy/treetime/treetime/treeanc.py#L880-L917)

Division by near-zero probabilities in the forward pass can amplify numerical errors for positions where a child contributes near-zero likelihood for some states. v0 avoids this by operating in log space (subtraction instead of division).

Golden master tests currently compare v1 dense marginal against v0 with tolerance 1e-6 to 1e-7, absorbing differences from the representation change.

A public Python report observed underflow while multiplying marginal subtree likelihoods for a large polytomy and state space [[issue](https://github.com/neherlab/treetime/issues/39)]. The maintainer identified log-space message accumulation as the remedy [[comment](https://github.com/neherlab/treetime/issues/39#issuecomment-383567696)]. The report is historical evidence for the numerical failure mode; the Rust forward-pass operations and reproduction conditions must be verified independently.
