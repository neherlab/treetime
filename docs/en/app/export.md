---
title: Export
description: Download output files and run packages, import a run package, and repeat a run with the command-line program.
---

## Download menu

"Download" in the run header lists, in this order:

- **Run package (.zip)**: the inputs, the settings, the log, and the outputs of the run; see [Run package](#run-package)
- **Main formats** of the analysis:
  - Time tree: Auspice JSON, the dated tree in Newick (branch lengths in years), and Nexus
  - Clock signal: the root-to-tip table and the rerooted tree in Newick
  - Ancestral sequences: the ancestral sequences in FASTA, the annotated tree in Nexus, and Auspice JSON
  - Optimize tree and Prune tree: the tree in Newick
- **All files...**: a dialog with every output file of the run and a description of each

The name of a downloaded file starts with the run name and the first characters of the run ID, so the files of different runs do not replace each other in your downloads folder.

When the browser can save to a chosen file (the page is served over HTTPS or from `localhost`, and the browser has the file save dialog of the File System Access API, as Chrome and Edge do), the web app asks where to save and writes the file to disk while it reads it. Otherwise the web app keeps the whole file in the memory of the page until the download starts, so a large run package can fail to download there.

## Run package

A run package is one zip file that holds a whole run: its record with the settings and the results summary, the log, the input files, the output files, the YAML configuration, and a shell script with the command line. Download it from the Download menu or from the menu of the run in the history.

A run package is the way to keep a run of the web app and to share a run: drop it on the start screen, or anywhere in the app, to import it. The imported run shows the same results, and you can make variants of it, as if you had run it in this app. The desktop app imports the packages of the web app, and the web app those of the desktop app. The app checks a package before the import: the record must be readable, the run must have ended, and each input file must have the size and the SHA-256 checksum that the record gives.

The package also repeats the run with the command-line program. Extract it and run the script in the extracted folder `treetime-<run-id>`:

```sh
unzip <package>.zip
cd treetime-<run-id>
sh command.sh
```

The script runs `treetime` with the configuration `run.yaml`, which reads the inputs from `inputs/` and writes the outputs into `out/`, with the seed of the run.

## Command line

"Command line" (keyboard shortcut C) shows the run or the open options as a `treetime` command or as a YAML configuration, with a copy button. The command holds the inputs, the settings that differ from the defaults, the outputs that the app adds to every run, and the seed of the run, with the paths of the run package. The dialog also gives:

- **Methods**: a paragraph for the methods section of a paper, with the analysis, the TreeTime version, and the key settings
- **Citation**: the paper to cite for the analysis
