# Smoke test fixtures

Inputs that `dev/smoke.toml` uses to reach failures and edge cases the example datasets do not contain. They are not examples: the smoke test does not run the `.yaml` file here as an example config, and no dataset directory lists these files.

| File                                         | Derived from           | Reaches                                                                                                                                                                                         |
| -------------------------------------------- | ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `zika-20-tree-without-branch-lengths.nwk`    | `zika/20/tree.nwk`     | Every branch length removed: `clock` on a tree without lengths                                                                                                                                  |
| `flu-h3n2-20-tree-single-child-polytomy.nwk` | `flu/h3n2/20/tree.nwk` | A single-child node `U` above the clade of `A/Oregon/15/2009` and `A/Hong_Kong/H090_695_V10/2009`, and the 15-tip clade collapsed into a polytomy: polytomy resolution with a single-child node |
| `clock-keep-root-and-reroot.yaml`            | -                      | A config that sets `keep_root` and `reroot`, which conflict on the command line                                                                                                                 |
