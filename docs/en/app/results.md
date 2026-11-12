---
title: Results
description: The run header, the overview, and the result views of each analysis.
---

A run page shows the result of one run. The options of the run stay in the panel at the left, so you can change an assumption next to the result; see [Settings of a shown run](/app/workspace#settings-of-a-shown-run).

## Run header

The header shows the run name, the status, the analysis, the dataset, the date of the run, and its duration. The pencil button next to the name renames the run. The actions of the header:

- **Download**: the run package and the main output files; see [Download menu](/app/export#download-menu)
- **Command line**: the command and the configuration that repeat the run; see [Command line](/app/export#command-line)
- **Create variant**: a new run from the settings of this run; see [Variants](/app/variants)
- **Menu**: rename, pin, compare with source, download run package, and delete

When the run wrote warnings, one line under the header gives their number and the first cause, such as "3 warnings: 3 tips dated outside their input date range". "Show" lists them all. A warning also shows next to the estimate that it affects.

## Views

The tabs under the header select the view. Each view exports its own data. The selection of a node and the filters of the tree stay when you change the view.

| View       | Content                                                                                                           | Analyses                                      |
| ---------- | ----------------------------------------------------------------------------------------------------------------- | --------------------------------------------- |
| Overview   | The main estimates and the tree                                                                                   | All                                           |
| Clock      | The root-to-tip regression, the residuals, the regression statistics, the clock outliers, and the undated samples | Time tree                                     |
| Population | The coalescent prior: the skyline or the coalescent time scale, and the coalescent settings                       | Time tree                                     |
| Model      | The substitution model and the entropy of each site                                                               | Time tree, Ancestral sequences, Optimize tree |
| Log        | The log of the run, by stage, with a filter, a search, and a copy action                                          | All                                           |
| Settings   | The settings of the run; "Changed only" shows only the settings that the run changed                              | All                                           |

A run that has not finished shows only Overview, Log, and Settings.

### Overview of a time tree

The estimates panel shows the root date with its interval and the clock rate with its standard deviation, such as "7.6e-4 per site per year, ± 4.1e-5", then the temporal signal and the number of samples used. "Run details" shows the clock model, the coalescent prior, the log likelihood, the number of iterations, and the duration.

When the run computed no interval for the dates, the panel says "No interval computed" and offers "Compute the interval": this starts a variant of the run with confidence intervals on.

The dated tree is the largest element. On a wide screen, the root-to-tip plot sits next to it. A click on a node in the tree or on a point in the plot selects the sample in both, and a side panel shows its details: its date, its clock data, its mutations, and why the clock fit excluded it, if it did.

### Overview of the other analyses

- **Clock signal**: the root-to-tip plot with the tree next to it, and the regression estimates
- **Ancestral sequences**: the tree colored by the mutations at the site that mutates on the most branches, the mutations of each branch, and the sites that mutated on several branches
- **Optimize tree** and **Prune tree**: the result tree, and a table of statistics of the input tree and the result tree

### Clock view

The Clock view shows the regression statistics as Tempest lists them: the clock rate, the root date from the x-intercept, the correlation coefficient r, R², and the residual mean square. The info button of each statistic gives its definition. The lists of the clock outliers and of the undated samples link to the tree and the plot.

## Tree viewer

A toolbar above the tree gives the coloring, the layout, the search, and the zoom to the root. "View options" opens all options of the Auspice viewer, among them the download of the tree as shown.

## Running and failed runs

While a run computes, the main area shows its stages, the clock rate of each iteration, and a cancel action. The options stay editable, so you can prepare the next run.

A failed run shows the error with its causes and the affected samples, and the stages of the log. "Edit and run again" opens the options at the setting that the error names.
