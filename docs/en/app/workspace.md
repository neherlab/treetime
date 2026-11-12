---
title: Workspace
description: The options panel, the input facts, the settings, the checks, and the preview before the first run.
---

The workspace shows the options of an analysis on the left and the main area on the right. The same layout serves a new analysis and a run: while you change an assumption, the result that motivated the change stays in view. Drag the border between the two panels to change their widths. On a narrow screen, the options open in a sheet, and the Run button stays in a bar at the bottom.

The header of the workspace names the dataset: the example name, or the file names of your data.

## Options panel

The options panel has four parts, from top to bottom.

### Analysis

Select the analysis: Time tree (selected by default), Clock signal, Ancestral sequences, Optimize tree, or Prune tree. The info button beside the selector describes the analysis and names its `treetime` command.

### Data

The input slots (tree, alignment, metadata) with their files, and the facts that TreeTime detected in them: the number of tips, the alignment length and alphabet, the metadata delimiter, the ID column, and the date column with the number of dates read. For example: "Dates: column `date`, 100 of 100 read as YYYY-MM-DD". "Change" beside a fact opens the setting that overrides the detection.

### Main settings

The key assumptions of the analysis, such as the clock rate, the root, the clock filter, the coalescent prior, and the confidence intervals. Each row shows the value that the run will use, also when you set nothing, such as "Automatic: least-squares rerooting". The help button of a row explains the setting and its options. The section title shows how many main settings differ from their defaults.

### All settings

"All settings" opens every setting of the analysis in a panel at the right, grouped by purpose:

- **Input data**: how TreeTime reads the inputs
- **Analysis**: the model and the inference, grouped by topic
- **Output format**: the settings that change the content of the output files, such as the annotations of the Newick output
- **Tree ordering and reproducibility**: the order of the nodes in the output tree, and the seed of the random steps

The search finds a setting by its name, its description, or its command-line flag. "Changed only" shows only the settings that differ from their defaults, and the reset button of a setting returns it to its default. The help of each setting gives its command-line flag, such as "Command line: `--clock-filter 4`". The app sets the output paths of each run itself, so they are not in the form.

### Checks and Run

The footer at the bottom of the panel checks the options while you edit them. It shows "Ready to run." or the findings: a finding that blocks the run, a warning, or advice. A finding names the settings by their labels and can offer a fix: a button that changes the settings, such as "Use covariation". Run starts the analysis; the keyboard shortcut is Ctrl+Enter (Cmd+Enter on macOS).

## Options menu

The menu button of the options panel offers:

- **Command line**: the options as a `treetime` command and as a YAML configuration; see [Command line](/app/export#command-line). The keyboard shortcut is C
- **Apply settings from a run**: the settings that a chosen run changed from the defaults replace the settings of the options. The inputs stay. The list shows the ended runs of the same analysis
- **Load config (YAML)**: load the settings of a YAML configuration file, such as one that `treetime` reads with `--config`

## Before the first run

Before you run, the main area shows what TreeTime can tell from the inputs:

- **Input data**: the facts of the tree, the alignment, and the metadata
- **Tree**: the input tree without dates, in divergence units
- **Root-to-tip regression**: the clock preview, the root-to-tip distance of each dated sample against its date, with the regression line and the clock-filter outliers in red. It uses the rooting and the clock-filter settings of the options, so you can check the temporal signal before the run

The preview computes trees of up to 5100 tips at once. For a larger tree, "Update preview" computes it on request.

## Settings of a shown run

A click on a run in the history shows its results and its options. The options are read-only until you change one: the first change turns them into a [variant](/app/variants) of the run. A bar above the results then lists the settings that differ from the shown run, such as "Settings changed since this run: coalescent prior, confidence. Run to update." Run starts the variant, and "Revert to this run's settings" puts back the options of the run.

"New analysis" in the sidebar starts from the default settings and no data. "Reset form" empties the form, and the Undo of its notice puts the form back.
