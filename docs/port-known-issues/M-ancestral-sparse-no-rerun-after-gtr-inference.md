# Sparse ancestral does not re-run marginal after GTR inference

In the ancestral command's sparse branch at [packages/treetime/src/commands/ancestral/run.rs#L138-L144](../../packages/treetime/src/commands/ancestral/run.rs#L138-L144), when `--model=infer` is used, the GTR is replaced after inference but the marginal reconstruction is not re-run with the new GTR. The dense branch at [run.rs#L193](../../packages/treetime/src/commands/ancestral/run.rs#L193) calls `update_marginal` a second time after GTR replacement.

## Impact

When `--model=infer` is used with sparse mode, the written `ancestral_sequences.fasta` contains sequences reconstructed using the initial dummy JC69 GTR, not the inferred GTR. The GTR stored on the partition is replaced, but the node profiles and reconstructed sequences still reflect the dummy model.

This is an unjustified divergence between the dense and sparse code paths. Dense mode correctly re-runs marginal after GTR inference; sparse mode does not.

## Affected code

- Sparse branch (missing re-run): [packages/treetime/src/commands/ancestral/run.rs#L138-L144](../../packages/treetime/src/commands/ancestral/run.rs#L138-L144)
- Dense branch (correct): [packages/treetime/src/commands/ancestral/run.rs#L193](../../packages/treetime/src/commands/ancestral/run.rs#L193) -- calls `update_marginal` after GTR replacement

## Fix

Add a second `update_marginal` call after GTR replacement in the sparse branch, matching the dense branch pattern. Verify that the resulting log-likelihood improves (or at least does not worsen) after the re-run.

## Related

- [M-core-dummy-gtr-initialization.md](M-core-dummy-gtr-initialization.md): the dummy GTR initialization pattern that creates this dependency
- [M-ancestral-dense-sparse-divergence.md](M-ancestral-dense-sparse-divergence.md): dense-sparse divergence from other sources
