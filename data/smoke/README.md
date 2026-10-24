# Smoke test fixtures

Inputs that `dev/smoke.toml` uses to reach failures and edge cases the example datasets do not contain. They are not examples: the smoke test does not run the `.yaml` file here as an example config, and no dataset directory lists these files.

| File                                               | Derived from                  | Reaches                                                                                                                                                                                         |
| -------------------------------------------------- | ----------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `zika-20-tree-without-internal-branch-lengths.nwk` | `zika/20/tree.nwk`            | Branch lengths of the internal edges removed: `clock` reads a missing length as 0                                                                                                               |
| `flu-h3n2-20-tree-single-child-polytomy.nwk`       | `flu/h3n2/20/tree.nwk`        | A single-child node `U` above the clade of `A/Oregon/15/2009` and `A/Hong_Kong/H090_695_V10/2009`, and the 15-tip clade collapsed into a polytomy: polytomy resolution with a single-child node |
| `mpox-clade-ii-1000-tree-optimized.nwk`            | `mpox/clade-ii/1000/tree.nwk` | Rerooted and optimized by an earlier `optimize` build, with `NODE_<n>` names: `optimize` and `timetree` on an input that an earlier run already rerooted                                        |
