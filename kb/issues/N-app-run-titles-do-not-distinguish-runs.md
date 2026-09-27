# Automatic run titles do not tell runs apart

The app names a run from its command and the flags that differ from the defaults. Runs of one command on different datasets with the same flags get the same title, and a long flag list fills the title. In the run list, the compare page and the "Date of this clade in other runs" chart, such runs show as the same truncated text, for example `Time tree, --output-selection 'nwk,nexu...`.

## Evidence

- `fn autoTitle()` in `packages/app-ui/src/settings/titles.ts` joins the command label with up to `TITLE_FLAGS` changed flags; the input files are not part of the title
- Six time-tree runs on the `flu/h3n2/20`, `dengue/100` and `mpox/clade-ii/1000` datasets carry the title `Time tree, --output-selection 'nwk,nexus,auspice,augur-node-data,gtr,reconstructed-nuc-fasta,clock-model,coalescent-tsv,tracelog,clock-csv'`, one of them with ` (edited)` appended

## Impact

- Users identify runs by date or by opening them, not by their title

## Possible changes

- Put the dataset (the tree or alignment file name) in the automatic title before the changed flags
- Show flag values only up to a length, so that the title keeps the flag names
