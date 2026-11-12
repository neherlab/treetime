---
title: Variants and comparison
description: Change one assumption of a run, read the difference to the source run, and compare two runs.
---

A variant is a run made from the settings of another run, its source, with one or more settings changed. A variant and its source form one line of work in the [run history](/app#run-history).

## Create a variant

"Create variant" in the run header makes the options of the run editable. The panel then names the source run and lists the settings that you changed, and the list grows while you edit. You can also start a variant with a change of any setting of a shown run; see [Settings of a shown run](/app/workspace#settings-of-a-shown-run). Run starts the variant. Its name lists the changed settings, such as "ebola/100, time tree: skyline".

A variant keeps the seed of its source run. TreeTime resolves polytomies with a random process, so two runs with different seeds can give slightly different estimates. With the same seed, a difference between the variant and its source comes from the changed settings only. To draw a new seed, clear the seed under All settings, in the section "Tree ordering and reproducibility".

## Difference to the source

The overview of a variant shows one line with the change of the main estimates against the source run, such as "Root date -46 d, clock rate -7.9 % against the source", with a link to the full comparison.

## Compare two runs

Open the comparison with "Compare with source" in the menu of a variant, or select two finished runs in the history and click "Compare". The swap button exchanges run A and run B.

- **Estimates**: the root date, the clock rate, the interval width, the log likelihood, the correlation coefficient r, and the excluded samples of both runs, with the difference in the column "Difference (B - A)". The comparison marks the log likelihoods as not comparable when the runs differ in their data, their coalescent prior, their substitution model, or the topology of their result trees, because then the two values belong to different models
- **Trees**: whether the trees have the same topology, and how far each shared ancestor moves in time between the two runs
- **Settings that differ**: the settings with different values in the two runs. The output settings that the app sets for each run are not in the list

When the two runs are different analyses, the page says so.
