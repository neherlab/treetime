---
title: The TreeTime app
description: Open a dataset, check its clock signal, run an analysis, and read the result in the web app or the desktop app.
---

The TreeTime app runs the analyses of the `treetime` command-line program on your own data or on an example dataset. Two apps share one interface:

- **Web app** ([treetime.neherlab.click](https://treetime.neherlab.click)): TreeTime runs inside the browser tab. Your files stay in the browser and nothing is uploaded to a server. The runs live in the tab; see [Runs in the web app and the desktop app](/app/runs)
- **Desktop app**: TreeTime runs on your computer and keeps every run in a runs folder on disk

Both apps compute with the same code as the command-line program, and every run gives its [command line](/app/export#command-line), so you can repeat it with `treetime`.

The usual work loop:

1. Open a dataset: drop your files or open an example
2. Look at the input facts and the clock preview, and correct what the app detected
3. Run the analysis and read the overview
4. Change one assumption in a [variant](/app/variants) and read the difference to the source run
5. [Export](/app/export) the result, or download the run package to keep the run

## Start screen

The start screen opens when you open the app at its address.

- **Drop zone**: drop the tree, the alignment, and the metadata here, or click "Choose files". The app reads the content of each file to find its kind: the file extension is only a hint. It puts each file into its input slot and opens the options of a new analysis. A file that the app cannot place stays out, and a notice gives the reason for each such file. Below the button, the app lists the formats that each input accepts
- **Configuration**: a YAML configuration in the drop zone loads its settings into the options, as "Load config (YAML)" in the [options menu](/app/workspace#options-menu) does
- **Run package**: a run package (`.zip`) in the drop zone imports the run with its results; see [Run package](/app/export#run-package)
- **Examples**: one card per example dataset, with the organism and the segment, the number of samples, the date range, and the sequence length. "Open" opens the dataset in a new analysis with the time tree selected. On the desktop app without examples, the gallery offers to download them
- **Recent work**: the latest lines of work, each with the name, the status, and the headline estimate of its newest run

You can also drop files anywhere else in the app: they go into the input slots of the open analysis.

## Run history

The sidebar lists the runs. A run and the variants made from it form one line of work: each variant shows under its source run. An entry shows the run name, a status icon, the headline estimate (such as the root date), the day of the run, and its size. "Filter by name" and the analysis filter narrow the list.

The menu of an entry has these actions:

- **Rename**: give the run your own name. The app names each run from the dataset, the analysis, and the settings that differ, such as "ebola/100, time tree: skyline, confidence"
- **Pin**: keep the line of work of the run at the top of the list
- **Compare with source**: open the [comparison](/app/variants#compare-two-runs) with the run that the variant came from
- **Download run package**: see [Run package](/app/export#run-package)
- **Delete**: see [Delete runs](/app/runs#delete-runs)

To act on several runs, choose "Select" in the menu of an entry, or Ctrl+click (Cmd+click on macOS) or Shift+click the entries. The footer of the sidebar then offers "Compare" for exactly two finished runs, "Delete", and "Clear".

"New analysis" at the top of the sidebar starts an analysis with the default settings and no data. The command palette (Ctrl+K, Cmd+K on macOS) finds runs, examples, and actions by name.
