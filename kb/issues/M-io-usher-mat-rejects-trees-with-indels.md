# UShER MAT output fails on every tree with an insertion or deletion

The UShER mutation-annotated tree writers (`mat-pb`, `mat-json`) stop the whole command when any branch carries an insertion or deletion:

```
Error:
   0: Node 'V20140729' has an insertion or deletion that UShER MAT cannot represent
```

`fn mat_mutation` rejects every mutation event that is not a substitution [packages/app-output/src/tree_output.rs#L376-L385](../../packages/app-output/src/tree_output.rs#L376-L385).

## Scope

`ancestral`, `optimize` and `timetree` fail with MAT output on 10 of the 15 small datasets: `dengue/20`, `dengue/100`, `ebola/20`, `ebola/100`, `lassa/L/20`, `lassa/L/50`, `mpox/clade-ii/20`, `mpox/clade-ii/100`, `rsv/a/20`, `rsv/a/100`. They succeed on `flu/h3n2/20`, `tb/20`, `tb/100`, `zika/20` and `zika/86`. `prune` fails on the same datasets with `--prune-empty` or `--merge-shared-mutations`; with `--prune-short` or `--prune-nodes-list` alone it writes MAT for `ebola/20`.

MAT output is not part of the default `--output-all` selection, so the failure appears only when MAT is selected, explicitly or through `--output-selection=all`. A user who asks for all outputs of a dataset with indels gets no output at all.

## Reference behavior

UShER has no gap state. `get_nuc_id` maps every character outside `ACGT` and the IUPAC ambiguity codes, the gap `-` included, to `0b1111`, the code of `N` (`Mutation_Annotated_Tree::get_nuc_id` in `src/mutation_annotated_tree.cpp` of [UShER](https://github.com/yatisht/usher), revision `ac9c982d`). UShER therefore stores a gap as missing data.

## Open question

How should the MAT writers handle insertions and deletions?

- Keep the error, and document that MAT output requires a tree without indels
- Store a deletion as missing data (`N`), as UShER does, and leave insertions out, with a warning that names the dropped events

## Smoke coverage

The `mat` rows of `dev/smoke.toml` run MAT output on every small dataset and declare the failing datasets. The default selections of `ancestral`, `optimize`, `prune` and `timetree` leave MAT out, so the other outputs of those datasets stay compared.
