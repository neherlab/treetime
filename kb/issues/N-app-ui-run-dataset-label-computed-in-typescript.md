# The run header computes the dataset label in TypeScript instead of reading it from Rust

Rust names the dataset of a run with `dataset_label()` (the example name, else the file name of the tree input), and builds the run name from it. The run header of the UI computes the same label again in TypeScript from the record, and the workspace header computes the label of a draft with its own rules. The apps take names and summaries from Rust, so a change of the Rust rule (for example a label for runs without a tree input) leaves the header with the old rule, and the run name and the header can name the dataset differently.

## Evidence

- Rust: `dataset_label()` [packages/app-backend/src/runs/run_name.rs#L38-L46](../../packages/app-backend/src/runs/run_name.rs#L38-L46); the record does not carry the label
- TypeScript copy for runs: `runDataset()` [packages/app-ui/src/runs/RunHeader.tsx#L64](../../packages/app-ui/src/runs/RunHeader.tsx#L64) (`record.example`, else the base name of the tree input)
- TypeScript rule for drafts: `datasetLabel()` [packages/app-ui/src/analysis/datasetLabel.ts#L6](../../packages/app-ui/src/analysis/datasetLabel.ts#L6) (the example folder, else the labels of all inputs)

## Fix direction

Add the dataset label to `RunRecord` or `RunSummary` from `dataset_label()`, and return the label of a draft from `check-config`, so both headers render a Rust value and `runDataset()` goes.

## Validation

- A route test checks the label of an example run and of a run with own files; the UI has no dataset rule of its own
