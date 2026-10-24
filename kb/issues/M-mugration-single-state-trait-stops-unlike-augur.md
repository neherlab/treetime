# Mugration stops on a trait with one state, while augur assigns that state

`treetime mugration` stops with an error when the trait column holds fewer than two known values:

```
A discrete trait needs at least 2 different states, found 1 state: denv2. Check that the trait column holds at least 2 different known values.
```

The check is in `fn StateAlphabet::new()` ([packages/treetime/src/alphabet/state_alphabet.rs#L15-L25](../../packages/treetime/src/alphabet/state_alphabet.rs#L15-L25)).

## Reference behavior

- **augur**: `augur traits` handles a single-state column before it calls TreeTime. It prints `WARNING: only one state found for discrete state reconstruction`, gives that state to every node, and writes no model (`augur/traits.py`, `def mugration_inference()`, the branch `len(unique_states)==1`). With no state, it gives every node no value
- **TreeTime 0.x**: `reconstruct_discrete_traits()` prints `mugration: only one or zero states found -- this doesn't make any sense` and returns no result (`packages/legacy/treetime/treetime/wrappers.py#L750-L752`)

v1 matches TreeTime 0.x. A Nextstrain build that runs `treetime mugration` in place of `augur traits` stops where augur goes on.

## Reproduction

`data/dengue/2000/mugration-serotype.yaml` reconstructs `serotype_genbank`, a trait of the Nextstrain dengue builds. Every dengue dataset in `data/` holds only DENV2, so the column has one state. `dev/smoke.toml` declares the expected failure of this config.

> [!IMPORTANT]
> **Decision required.** Choose the behavior for a trait with one known state:
>
> - Give the state to every node with probability 1, write the trait outputs, and write no model, as augur does. The run warns that the trait has one state
> - Keep the error, as TreeTime 0.x stops too. Pipelines that replace `augur traits` then leave the step out
>
> A trait with no known value is a separate case: augur gives every node no value, and the error stays useful there.
