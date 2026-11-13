# Firefox and Safari hold a whole web2 download in page memory

`saveRunOutput()` streams a run package or an output file to disk through `showSaveFilePicker()`. Firefox and Safari have no `showSaveFilePicker()`, so there the page collects the whole download as one Blob before the save starts. A run package can be larger than 1 GB (ancestral on `mpox/clade-ii/2000` with default outputs is about 1.18 GB), so the download needs about twice its size in page memory, next to the memory of the back-end worker that holds the runs.

## Evidence

- `saveRunOutput()` and `browserFilePicker()` in `packages/app-ui/src/runs/saveRunOutput.ts`: the stream path when the picker exists, the Blob path otherwise
- Measured in Bun by driving `saveRunOutput()` with a body of 1126 MiB: the Blob path raised the peak memory by 2256 MiB, the stream path by 20 MiB. Browser memory was not measured

> [!IMPORTANT]
> **Investigation required.** Measure the download of a 1 GB run package in Firefox and Safari, and find the size at which it fails. Options afterwards: a service worker that streams the download to the browser's download manager, or a size note in the Download menu of these browsers.

## Validation

- A run package of the largest example downloads in Firefox and Safari, or the app tells the user before the download that it cannot
