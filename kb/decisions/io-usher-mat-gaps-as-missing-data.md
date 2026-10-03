# UShER MAT output writes gaps as missing data

The UShER mutation-annotated tree writers (`mat-pb`, `mat-json`) write a deletion as missing data (`N`) and leave out insertions that the MAT reference cannot represent, with one warning per run that counts what was written as missing data or left out. A tree with insertions or deletions therefore no longer stops the whole command.

**Type**: Behavior change in v1 (no v0 counterpart: v0 writes no MAT).

**Status**: Approved by the maintainer team. The rejected alternative kept the error for every tree with an insertion or deletion and documented that MAT output requires a tree without indels. That error stopped `ancestral`, `optimize`, `prune` and `timetree` with MAT output on 10 of the 15 small datasets (`dengue`, `ebola`, `lassa`, `mpox`, `rsv`), so a user who selected all outputs got none.

**v1 location**: `fn mat_from_graph()` [packages/app-output/src/tree_output.rs#L407](../../packages/app-output/src/tree_output.rs#L407) converts the indels with `fn gaps_as_missing_data()` [packages/app-output/src/tree_output.rs#L472](../../packages/app-output/src/tree_output.rs#L472), passes them through the missing-data bridge `UnknownBridge` [packages/app-output/src/mutation_filter.rs#L108](../../packages/app-output/src/mutation_filter.rs#L108), and counts them in `MatGapCounts` [packages/app-output/src/tree_output.rs#L356](../../packages/app-output/src/tree_output.rs#L356). `fn converted_mat()` [packages/app-output/src/tree_output.rs#L109](../../packages/app-output/src/tree_output.rs#L109) converts the tree once for both MAT formats and emits the warning.

## Reference behavior

UShER has no gap state:

- `Mutation_Annotated_Tree::get_nuc_id` maps every character outside `ACGT` and the IUPAC ambiguity codes, the gap `-` included, to `0b1111`, the code of `N` (`src/mutation_annotated_tree.cpp` of [UShER](https://github.com/yatisht/usher), revision `ac9c982d`). UShER therefore reads a gap as missing data
- The protobuf reader decodes `ref_nuc` and `par_nuc` as one base each (`1 << ref_nuc`, `1 << par_nuc`, same file). A column in which the reference has no base cannot carry a mutation
- The reader prints `WARNING: Mutations not sorted!` and sorts a node's mutations when they are not in position order (same file)

## Rule

The MAT reference is the root sequence. For each edge, on the nucleotide track:

- **Deletion**: every deleted site becomes a substitution from the parent state to `N`. The MAT writer already stores `N` as missing data: `UnknownBridge` hides the change and remembers the state before it
- **Insertion below a deletion**: every inserted site becomes a substitution from `N` to the inserted state. `UnknownBridge` reports it against the state before the deletion, or reports nothing when the inserted state equals that state
- **Deleted unknown site**: a deleted `N` hides no state, so a later insertion at that site has nothing to restore and is not reported
- **Root gap column**: a column in which the root sequence has a gap has no reference base. Every insertion and substitution in that column is left out
- **Order**: the mutations of each node are sorted by position
- **Warning**: when any event was converted or left out, the command logs one warning per run, for example `UShER MAT has no gap state: wrote 25 deletion(s) as missing data (N)`. It counts deletion events, insertion events with a site in a root gap column, and substitutions in root gap columns

Amino-acid mutations are outside this rule; MAT still rejects them ([kb/issues/M-io-usher-mat-rejects-amino-acid-mutations.md](../issues/M-io-usher-mat-rejects-amino-acid-mutations.md)).

## Example

Tree `((C,D)B,E)root`:

| Root sequence | Edge events                                    | MAT mutations | Warning counts                       |
| ------------- | ---------------------------------------------- | ------------- | ------------------------------------ |
| `ACGT`        | B deletes `CG` at columns 2-3                  | none          | 1 deletion                           |
| `ACGT`        | B deletes `C` at column 2, C inserts `T` there | C: `C2T`      | 1 deletion                           |
| `A-GT`        | B inserts `C` at column 2, C substitutes `C2T` | none          | 1 insertion, 1 substitution left out |

The tests `test_tree_output_mat_writes_gaps_as_missing_data` ([packages/app-output/src/**tests**/test_tree_output_mat.rs](../../packages/app-output/src/__tests__/test_tree_output_mat.rs)) check these cases. `test_mat_gaps_written_as_missing_data_with_one_warning` ([packages/app-commands/src/commands/ancestral/**tests**/test_mat_gaps.rs](../../packages/app-commands/src/commands/ancestral/__tests__/test_mat_gaps.rs)) runs `ancestral` in dense and sparse mode on an alignment with a deletion and a root-gap insertion, and checks the MAT mutations and the single warning.

## Why

- UShER itself reads a gap as missing data, so a MAT that TreeTime writes reads back in UShER with the same meaning
- The MAT writer already handles `N` with `UnknownBridge`. A deletion written as `N` goes through the same mechanism, so a base that is deleted and inserted again keeps a correct parent state
- One output format that cannot represent indels no longer prevents every other output of the run

## Impact

- MAT output is now written for the datasets with indels. `ancestral` on `rsv/a/20` warns `wrote 25 deletion(s) as missing data (N)`
- MAT output loses the indels. Newick comments, Auspice JSON, augur node-data JSON and the reconstructed FASTA keep them
- The mutations of a node are now written in position order on every dataset
