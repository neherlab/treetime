# Application UI ignores inputs and displays synthetic scientific results

Visible parameter controls keep only component-local state and do not affect requests. Dataset switching retains absent files, command results are discarded, fixed mock scientific values are displayed after completion, and caught errors are hidden. The UI records only the name and size of a browser-selected file, while the server expects such files to be uploaded into a run's `inputs/` folder before the run starts.

## Evidence

- `function ParamForm()` keeps controls in component state without propagating them into command requests [packages/app-ui/src/components/ParamForm.tsx#L308-L364](../../packages/app-ui/src/components/ParamForm.tsx#L308-L364).
- `function RunButton()` constructs requests from fixed fields and discards successful result values [packages/app-ui/src/components/RunButton.tsx#L19-L76](../../packages/app-ui/src/components/RunButton.tsx#L19-L76).
- `function ResultsPanel()` renders fixed clock and optimization values [packages/app-ui/src/components/ResultsPanel.tsx#L89-L97](../../packages/app-ui/src/components/ResultsPanel.tsx#L89-L97).

## Impact

The UI can run a different computation from the one configured by the user and then present fabricated values as its result. This violates scientific-output integrity even when the backend command itself is correct.

## Potential solutions

- O1. Connect typed UI state, real results, and visible errors end to end while accepting only server-discovered dataset paths.
- O2. Keep the execution UI read-only until those contracts exist and remove controls/results that imply unsupported behavior.

## Recommendation

Store typed command state centrally, replace dataset selections atomically, preserve typed command results, render only real result data, and show causal errors. Until a real result projection exists, keep that result view unavailable after execution. Browser-selected files are uploaded into the run's `inputs/` folder, and the run configuration names the uploaded copies.
