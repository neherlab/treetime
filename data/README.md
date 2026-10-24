# Example configs and pipelines

This directory holds the example datasets and, next to each one, ready-to-run configuration files for TreeTime. Every command reads its full argument object from a file with `--config`, and the `pipeline` command runs an ordered chain of commands from a single file. The examples double as living documentation: each file is a saved, reviewed, version-controlled invocation you can copy and adapt.

Examples run from any folder. A relative path inside a file resolves from the folder of that file, so a dataset folder works wherever it is unpacked. A relative path given as a command-line flag resolves from the working directory. The files name no output paths: choose the output folder with `--output-all`, for example `tmp/` (git-ignored) when running from the repository root.

## Per-command configs

A per-command config is exactly the command's argument object. Its keys map one-to-one to the CLI flags: `--alignment a.fasta` is `alignment: [a.fasta]`, `--clock-rate 0.003` is `clock_rate: 0.003`. An explicit flag on the command line still overrides the file (flag > config > default).

| File | Command | Demonstrates |
| --- | --- | --- |
| [`dengue/100/clock.yaml`](dengue/100/clock.yaml) | `clock` | Root-to-tip clock regression |
| [`dengue/2000/ancestral.yaml`](dengue/2000/ancestral.yaml) | `ancestral` | Ancestral sequences with the default settings |
| [`dengue/2000/clock.yaml`](dengue/2000/clock.yaml) | `clock` | Root-to-tip regression with the default settings |
| [`dengue/2000/homoplasy.yaml`](dengue/2000/homoplasy.yaml) | `homoplasy` | Homoplasies with the default settings |
| [`dengue/2000/mugration-country.yaml`](dengue/2000/mugration-country.yaml) | `mugration` | Ancestral sampling country |
| [`dengue/2000/mugration-serotype.yaml`](dengue/2000/mugration-serotype.yaml) | `mugration` | Ancestral serotype of 1954 dengue genomes |
| [`dengue/2000/optimize-dense.yaml`](dengue/2000/optimize-dense.yaml) | `optimize` | Branch-length optimization in the dense representation |
| [`dengue/2000/optimize.yaml`](dengue/2000/optimize.yaml) | `optimize` | Branch-length optimization with the default settings |
| [`dengue/2000/prune-collapse.yaml`](dengue/2000/prune-collapse.yaml) | `prune` | Tree cleanup: short, empty, and shared-mutation branches |
| [`dengue/2000/timetree.yaml`](dengue/2000/timetree.yaml) | `timetree` | Time tree with the default settings |
| [`ebola/20/ancestral-parsimony.yaml`](ebola/20/ancestral-parsimony.yaml) | `ancestral` | Ancestral reconstruction by Fitch parsimony |
| [`ebola/20/prune-remove.yaml`](ebola/20/prune-remove.yaml) | `prune` | Remove named tips from the tree |
| [`ebola/362/ancestral.yaml`](ebola/362/ancestral.yaml) | `ancestral` | Ancestral sequences with the default settings |
| [`ebola/362/clock.yaml`](ebola/362/clock.yaml) | `clock` | Root-to-tip regression with the default settings |
| [`ebola/362/homoplasy.yaml`](ebola/362/homoplasy.yaml) | `homoplasy` | Homoplasies with the default settings |
| [`ebola/362/mugration-country.yaml`](ebola/362/mugration-country.yaml) | `mugration` | Ancestral country of the 2014-2016 West African Ebola epidemic |
| [`ebola/362/mugration-division.yaml`](ebola/362/mugration-division.yaml) | `mugration` | Ancestral district of the 2014-2016 West African Ebola epidemic |
| [`ebola/362/optimize.yaml`](ebola/362/optimize.yaml) | `optimize` | Branch-length optimization with the default settings |
| [`ebola/362/prune-collapse.yaml`](ebola/362/prune-collapse.yaml) | `prune` | Tree cleanup: short, empty, and shared-mutation branches |
| [`ebola/362/timetree-dense.yaml`](ebola/362/timetree-dense.yaml) | `timetree` | Time tree with a skyline coalescent in the dense representation |
| [`ebola/362/timetree.yaml`](ebola/362/timetree.yaml) | `timetree` | Time tree with the default settings |
| [`flu/h3n2/20/optimize-reroot-min-dev.yaml`](flu/h3n2/20/optimize-reroot-min-dev.yaml) | `optimize` | Branch-length optimization with divergence-based rerooting |
| [`flu/h3n2/200/timetree-fixed-clock.yaml`](flu/h3n2/200/timetree-fixed-clock.yaml) | `timetree` | Time-scaled tree with a fixed clock rate and date confidence intervals |
| [`flu/h3n2/500/ancestral-codon-partitions.yaml`](flu/h3n2/500/ancestral-codon-partitions.yaml) | `ancestral` | Ancestral sequences with one substitution model per codon position |
| [`flu/h3n2/500/ancestral.yaml`](flu/h3n2/500/ancestral.yaml) | `ancestral` | Ancestral sequences with the default settings |
| [`flu/h3n2/500/clock.yaml`](flu/h3n2/500/clock.yaml) | `clock` | Root-to-tip regression with the default settings |
| [`flu/h3n2/500/homoplasy.yaml`](flu/h3n2/500/homoplasy.yaml) | `homoplasy` | Homoplasies with the default settings |
| [`flu/h3n2/500/mugration-location.yaml`](flu/h3n2/500/mugration-location.yaml) | `mugration` | Ancestral sampling country |
| [`flu/h3n2/500/mugration-region.yaml`](flu/h3n2/500/mugration-region.yaml) | `mugration` | Ancestral sampling region |
| [`flu/h3n2/500/optimize-codon-partitions.yaml`](flu/h3n2/500/optimize-codon-partitions.yaml) | `optimize` | Branch lengths with one substitution model per codon position |
| [`flu/h3n2/500/optimize.yaml`](flu/h3n2/500/optimize.yaml) | `optimize` | Branch-length optimization with the default settings |
| [`flu/h3n2/500/prune-collapse.yaml`](flu/h3n2/500/prune-collapse.yaml) | `prune` | Tree cleanup: short, empty, and shared-mutation branches |
| [`flu/h3n2/500/timetree-codon-partitions.yaml`](flu/h3n2/500/timetree-codon-partitions.yaml) | `timetree` | Time tree with one substitution model per codon position |
| [`flu/h3n2/500/timetree-skyline-confidence.yaml`](flu/h3n2/500/timetree-skyline-confidence.yaml) | `timetree` | Influenza H3N2 time tree with an estimated clock and a skyline coalescent |
| [`flu/h3n2/500/timetree.yaml`](flu/h3n2/500/timetree.yaml) | `timetree` | Time tree with the default settings |
| [`lassa/L/20/timetree.yaml`](lassa/L/20/timetree.yaml) | `timetree` | Lassa virus time tree with the default settings |
| [`lassa/L/500/ancestral.yaml`](lassa/L/500/ancestral.yaml) | `ancestral` | Ancestral sequences with the default settings |
| [`lassa/L/500/clock.yaml`](lassa/L/500/clock.yaml) | `clock` | Root-to-tip regression with the default settings |
| [`lassa/L/500/homoplasy.yaml`](lassa/L/500/homoplasy.yaml) | `homoplasy` | Homoplasies with the default settings |
| [`lassa/L/500/mugration-country.yaml`](lassa/L/500/mugration-country.yaml) | `mugration` | Ancestral sampling country |
| [`lassa/L/500/optimize.yaml`](lassa/L/500/optimize.yaml) | `optimize` | Branch-length optimization with the default settings |
| [`lassa/L/500/prune-collapse.yaml`](lassa/L/500/prune-collapse.yaml) | `prune` | Tree cleanup: short, empty, and shared-mutation branches |
| [`lassa/L/500/timetree-relaxed-clock.yaml`](lassa/L/500/timetree-relaxed-clock.yaml) | `timetree` | Lassa virus time tree with a relaxed clock |
| [`lassa/L/500/timetree.yaml`](lassa/L/500/timetree.yaml) | `timetree` | Time tree with the default settings |
| [`mpox/clade-ii/2000/ancestral.yaml`](mpox/clade-ii/2000/ancestral.yaml) | `ancestral` | Ancestral sequences with the default settings |
| [`mpox/clade-ii/2000/clock-outgroup.yaml`](mpox/clade-ii/2000/clock-outgroup.yaml) | `clock` | Root-to-tip regression on the outgroup root of the Nextstrain hmpxv1 build |
| [`mpox/clade-ii/2000/clock.yaml`](mpox/clade-ii/2000/clock.yaml) | `clock` | Root-to-tip regression with the default settings |
| [`mpox/clade-ii/2000/homoplasy.yaml`](mpox/clade-ii/2000/homoplasy.yaml) | `homoplasy` | Homoplasies with the default settings |
| [`mpox/clade-ii/2000/mugration-country.yaml`](mpox/clade-ii/2000/mugration-country.yaml) | `mugration` | Ancestral sampling country |
| [`mpox/clade-ii/2000/optimize.yaml`](mpox/clade-ii/2000/optimize.yaml) | `optimize` | Branch-length optimization with the default settings |
| [`mpox/clade-ii/2000/prune-collapse.yaml`](mpox/clade-ii/2000/prune-collapse.yaml) | `prune` | Tree cleanup: short, empty, and shared-mutation branches |
| [`mpox/clade-ii/2000/timetree.yaml`](mpox/clade-ii/2000/timetree.yaml) | `timetree` | Time tree with the default settings |
| [`rsv/a/100/prune-collapse.yaml`](rsv/a/100/prune-collapse.yaml) | `prune` | Tree cleanup: drop short and uninformative branches, collapse shared mutations |
| [`rsv/a/2000/ancestral-dense.yaml`](rsv/a/2000/ancestral-dense.yaml) | `ancestral` | Marginal reconstruction in the dense representation |
| [`rsv/a/2000/ancestral.yaml`](rsv/a/2000/ancestral.yaml) | `ancestral` | Ancestral sequences with the default settings |
| [`rsv/a/2000/clock.yaml`](rsv/a/2000/clock.yaml) | `clock` | Root-to-tip regression with the default settings |
| [`rsv/a/2000/homoplasy.yaml`](rsv/a/2000/homoplasy.yaml) | `homoplasy` | Homoplasies with the default settings |
| [`rsv/a/2000/mugration-country.yaml`](rsv/a/2000/mugration-country.yaml) | `mugration` | Ancestral sampling country |
| [`rsv/a/2000/optimize.yaml`](rsv/a/2000/optimize.yaml) | `optimize` | Branch-length optimization with the default settings |
| [`rsv/a/2000/prune-collapse.yaml`](rsv/a/2000/prune-collapse.yaml) | `prune` | Tree cleanup: short, empty, and shared-mutation branches |
| [`rsv/a/2000/timetree-skyline-confidence.yaml`](rsv/a/2000/timetree-skyline-confidence.yaml) | `timetree` | RSV-A time tree with a skyline coalescent and date intervals |
| [`rsv/a/2000/timetree.yaml`](rsv/a/2000/timetree.yaml) | `timetree` | Time tree with the default settings |
| [`sc2/2844/ancestral.yaml`](sc2/2844/ancestral.yaml) | `ancestral` | Marginal ancestral reconstruction on 2504 SARS-CoV-2 genomes |
| [`sc2/5100/ancestral.yaml`](sc2/5100/ancestral.yaml) | `ancestral` | Ancestral sequences with the default settings |
| [`sc2/5100/clock-outgroup.yaml`](sc2/5100/clock-outgroup.yaml) | `clock` | Root-to-tip regression on the Wuhan-Hu-1 root |
| [`sc2/5100/clock.yaml`](sc2/5100/clock.yaml) | `clock` | Root-to-tip regression with the default settings |
| [`sc2/5100/homoplasy.yaml`](sc2/5100/homoplasy.yaml) | `homoplasy` | Homoplasies with the default settings |
| [`sc2/5100/mugration-country.yaml`](sc2/5100/mugration-country.yaml) | `mugration` | Ancestral geography over 95 countries on 5100 SARS-CoV-2 genomes |
| [`sc2/5100/mugration-region.yaml`](sc2/5100/mugration-region.yaml) | `mugration` | Ancestral region of 5100 SARS-CoV-2 genomes |
| [`sc2/5100/optimize.yaml`](sc2/5100/optimize.yaml) | `optimize` | Branch-length optimization with the default settings |
| [`sc2/5100/prune-collapse.yaml`](sc2/5100/prune-collapse.yaml) | `prune` | Tree cleanup: short, empty, and shared-mutation branches |
| [`sc2/5100/timetree-input-branch-lengths.yaml`](sc2/5100/timetree-input-branch-lengths.yaml) | `timetree` | Time tree from the input branch lengths, without sequences |
| [`sc2/5100/timetree-skyline-confidence.yaml`](sc2/5100/timetree-skyline-confidence.yaml) | `timetree` | SARS-CoV-2 time tree with a skyline coalescent, 5100 genomes |
| [`sc2/5100/timetree.yaml`](sc2/5100/timetree.yaml) | `timetree` | Time tree with the default settings |
| [`tb/149/ancestral-vcf.yaml`](tb/149/ancestral-vcf.yaml) | `ancestral` | Ancestral reconstruction from the VCF of the Lee 2015 TB isolates |
| [`tb/149/ancestral.yaml`](tb/149/ancestral.yaml) | `ancestral` | Ancestral sequences with the default settings |
| [`tb/149/clock.yaml`](tb/149/clock.yaml) | `clock` | Root-to-tip regression with the default settings |
| [`tb/149/homoplasy-constant-sites.yaml`](tb/149/homoplasy-constant-sites.yaml) | `homoplasy` | Homoplasies of TB isolates from an informative-sites alignment |
| [`tb/149/homoplasy.yaml`](tb/149/homoplasy.yaml) | `homoplasy` | Homoplasies with the default settings |
| [`tb/149/mugration-cluster.yaml`](tb/149/mugration-cluster.yaml) | `mugration` | Ancestral cluster of Mycobacterium tuberculosis isolates |
| [`tb/149/optimize.yaml`](tb/149/optimize.yaml) | `optimize` | Branch-length optimization with the default settings |
| [`tb/149/prune-collapse.yaml`](tb/149/prune-collapse.yaml) | `prune` | Tree cleanup: short, empty, and shared-mutation branches |
| [`tb/149/timetree.yaml`](tb/149/timetree.yaml) | `timetree` | Time tree with the default settings |
| [`zika/20/ancestral.yaml`](zika/20/ancestral.yaml) | `ancestral` | Marginal ancestral reconstruction with an inferred GTR model |
| [`zika/20/mugration-country-weights.yaml`](zika/20/mugration-country-weights.yaml) | `mugration` | Discrete ancestral geography ("mugration") over sampling country |
| [`zika/86/ancestral-phylip.yaml`](zika/86/ancestral-phylip.yaml) | `ancestral` | Ancestral reconstruction from a PHYLIP alignment |
| [`zika/86/ancestral.yaml`](zika/86/ancestral.yaml) | `ancestral` | Ancestral sequences with the default settings |
| [`zika/86/clock.yaml`](zika/86/clock.yaml) | `clock` | Root-to-tip regression with the default settings |
| [`zika/86/homoplasy.yaml`](zika/86/homoplasy.yaml) | `homoplasy` | Homoplasies with the default settings |
| [`zika/86/mugration-country.yaml`](zika/86/mugration-country.yaml) | `mugration` | Ancestral sampling country |
| [`zika/86/mugration-region.yaml`](zika/86/mugration-region.yaml) | `mugration` | Ancestral sampling region |
| [`zika/86/optimize.yaml`](zika/86/optimize.yaml) | `optimize` | Branch-length optimization with the default settings |
| [`zika/86/prune-collapse.yaml`](zika/86/prune-collapse.yaml) | `prune` | Tree cleanup: short, empty, and shared-mutation branches |
| [`zika/86/timetree.yaml`](zika/86/timetree.yaml) | `timetree` | Time tree with the default settings |

Run one with, for example:

```bash
treetime ancestral --config data/zika/20/ancestral.yaml --output-all tmp/zika-ancestral
```

## Pipelines

A pipeline file lists named steps, each an analysis command with its arguments. Steps run sequentially in one process. A step reads an earlier step's output with `{{ steps.<name>.outputs.<selection> }}`, so intermediate paths are never hand-written. `vars` names shared values once. `treetime pipeline --output-all <dir>`, or a top-level `output_all` in the file, gives each step its own `<dir>/<step>/` output directory; the flag takes precedence.

| File | Chain | Demonstrates |
| --- | --- | --- |
| [`dengue/2000/pipeline-full.yaml`](dengue/2000/pipeline-full.yaml) | `optimize` -> `ancestral` -> `timetree` -> `mugration` | Dengue build with a skyline coalescent and ancestral geography |
| [`dengue/2000/pipeline-nextstrain.yaml`](dengue/2000/pipeline-nextstrain.yaml) | `timetree` -> `ancestral` -> `mugration` -> `mugration` | Nextstrain dengue denv2 genome build, TreeTime part |
| [`ebola/362/pipeline-nextstrain.yaml`](ebola/362/pipeline-nextstrain.yaml) | `timetree` -> `ancestral` -> `mugration` -> `mugration` | Nextstrain ebola west-africa-2014 build, TreeTime part |
| [`flu/h3n2/200/pipeline-fixed-clock.yaml`](flu/h3n2/200/pipeline-fixed-clock.yaml) | `optimize` -> `timetree` | Optimize then date, with a fixed clock rate defined once |
| [`flu/h3n2/500/pipeline-nextstrain.yaml`](flu/h3n2/500/pipeline-nextstrain.yaml) | `optimize` -> `prune` -> `timetree` -> `ancestral` -> `mugration` | Nextstrain seasonal-flu h3n2 NA build, TreeTime part |
| [`lassa/L/500/pipeline-nextstrain.yaml`](lassa/L/500/pipeline-nextstrain.yaml) | `timetree` -> `ancestral` -> `mugration` | Nextstrain lassa L segment build, TreeTime part |
| [`mpox/clade-ii/2000/pipeline-full.yaml`](mpox/clade-ii/2000/pipeline-full.yaml) | `optimize` -> `ancestral` -> `timetree` -> `mugration` | All four analysis commands on 2000 clade II mpox genomes |
| [`mpox/clade-ii/2000/pipeline-nextstrain-new-tt.yaml`](mpox/clade-ii/2000/pipeline-nextstrain-new-tt.yaml) | `optimize` -> `ancestral` -> `timetree` | Nextstrain mpox WIP/new_tt build with TreeTime v1 |
| [`mpox/clade-ii/2000/pipeline-nextstrain.yaml`](mpox/clade-ii/2000/pipeline-nextstrain.yaml) | `timetree` -> `ancestral` | Nextstrain mpox hmpxv1 build, TreeTime part |
| [`rsv/a/2000/pipeline-nextstrain.yaml`](rsv/a/2000/pipeline-nextstrain.yaml) | `timetree` -> `ancestral` -> `mugration` -> `mugration` | Nextstrain rsv a genome build, TreeTime part |
| [`sc2/2844/pipeline-undated.yaml`](sc2/2844/pipeline-undated.yaml) | `optimize` -> `ancestral` | SARS-CoV-2 branch-length optimization and ancestral reconstruction |
| [`sc2/5100/pipeline-nextstrain.yaml`](sc2/5100/pipeline-nextstrain.yaml) | `timetree` -> `ancestral` -> `mugration` | Nextstrain ncov global build, TreeTime part |
| [`zika/20/pipeline-dated-geography.yaml`](zika/20/pipeline-dated-geography.yaml) | `timetree` -> `mugration` | Phylodynamics + phylogeography in one file: time tree, then ancestral geography |
| [`zika/86/pipeline-geography-map.yaml`](zika/86/pipeline-geography-map.yaml) | `timetree` -> `mugration` | Zika virus spread in the Americas: time tree, then weighted ancestral geography |
| [`zika/86/pipeline-nextstrain.yaml`](zika/86/pipeline-nextstrain.yaml) | `timetree` -> `ancestral` -> `mugration` -> `mugration` | Nextstrain zika build, TreeTime part |

Run a pipeline, preview its plan, or run only part of it:

```bash
treetime pipeline --config data/zika/20/pipeline-dated-geography.yaml --output-all tmp/zika-pipeline
treetime pipeline --config data/zika/20/pipeline-dated-geography.yaml --output-all tmp/zika-pipeline --dry-run
treetime pipeline --config data/zika/20/pipeline-dated-geography.yaml --output-all tmp/zika-pipeline --steps mugration
```

The mugration step of `zika/86/pipeline-geography-map.yaml` reads two tables of the sampling countries. [`zika/country_weights.csv`](zika/country_weights.csv) gives the equilibrium probability of each country. [`zika/country_lat_longs.tsv`](zika/country_lat_longs.tsv) gives its coordinates in the format of augur's `--lat-longs`, so the Auspice tree of the step shows the countries on a map. The coordinates are the label points of the countries in Natural Earth 1:50m Admin 0 Map Units, version 5.1.2 (public domain), which list French Guiana, Guadeloupe and Martinique as units of their own.

`--dry-run` resolves the whole plan (each step's inputs, its concrete output files, and every chained path) and prints it without running anything. `--steps` runs a named subset in list order; a referenced upstream step that is not selected must already have its outputs on disk. When a step fails, the error names the completed steps and prints the exact `--steps` value that resumes the rest.

## Names and places

Every config lives in the dataset folder whose files it reads, and a file that several sizes of a dataset share lives one folder up, as `zika/country_weights.csv` does. The largest dataset of each family has a config for every command; the small datasets keep the quick examples that `just smoke` runs.

- `<command>.yaml` runs the defaults of the command and sets only what the dataset needs, such as its id column. A command that needs a choice names it: `mugration-<attribute>.yaml` and `prune-<operation>.yaml`
- `<command>-<variant>.yaml` names the setting that makes the variant: `ancestral-parsimony.yaml`, `optimize-dense.yaml`, `timetree-relaxed-clock.yaml`, `timetree-skyline-confidence.yaml`
- `pipeline-nextstrain.yaml` runs the TreeTime part of a Nextstrain build of the pathogen, with the workflow commit and its augur settings in the header; `pipeline-full.yaml` runs all four analysis commands; other pipelines name their purpose, such as `pipeline-dated-geography.yaml`
- A config that fails today stays in place, and `dev/smoke.toml` names the known issue under `[config.fail]`. The header of the file names the issue too
- A per-command config names each input file directly, because the apps upload the files that the index lists for it. Protein reconstruction from per-CDS translations, whose setting is a file name template, therefore appears only in the Nextstrain pipelines of mpox and RSV-A

## Why a config or pipeline instead of a shell command

- One file, not one long flag string. The invocation is reviewable in a pull request, diffable across runs, and stored next to the data it analyzes.
- No repetition. `vars` holds a shared value such as the clock rate once; every step references it. A whole-value reference such as `{{ vars.clock_rate }}` keeps its JSON type, so a number stays a number.
- No path bookkeeping. `{{ steps.optimize.outputs.nwk }}` resolves to the optimized tree the first step actually writes. Rename an output directory in one place and the chain follows.
- Errors before computation. `--dry-run` resolves and type-checks the whole plan up front. A misspelled key is rejected with a suggestion (`unknown field opt_methd`; `unknown top-level key step; did you mean steps?`) rather than running silently or halfway.
- One process. Steps share the thread pool; there is no shell loop, no temporary files to name, and no second tool to install.

## Why a pipeline instead of Snakemake for the TreeTime part

- No workflow engine to stand up. A Snakemake build needs a `Snakefile`, a `config.yaml`, wildcards, and one rule per tool. The pipeline file needs only the TreeTime binary.
- Chaining by output identity, not by string. A Snakemake rule wires steps with literal `input:`/`output:` path strings that must match each tool's real filenames. The pipeline references an output by its selection tag (`nwk`, `auspice`, ...) and resolves the concrete path itself.
- Validation is static and typed. `--dry-run` and the JSON-schema-backed argument objects catch a bad parameter or a dangling step reference before anything runs; a Snakemake rule only surfaces a bad TreeTime parameter when that rule executes.
- Reproducible and portable. The config is the exact serialized argument shape, with one dependency, committed beside the data.

This is not a replacement for a full multi-sample workflow manager: incremental rebuilds, cluster scheduling, and non-TreeTime rules (subsampling, alignment, tree building, export) still belong to Snakemake and Augur. The pipeline replaces the TreeTime glue inside such a build. On the Nextstrain mpox `WIP/new_tt` phylogenetic workflow, for instance, Snakemake rules run `treetime optimize` in place of `augur refine` and `treetime ancestral` in place of `augur ancestral` and `augur translate`. [`mpox/clade-ii/2000/pipeline-nextstrain-new-tt.yaml`](mpox/clade-ii/2000/pipeline-nextstrain-new-tt.yaml) runs these steps and a dating step in one file, and every `pipeline-nextstrain.yaml` runs the TreeTime part of a Nextstrain build the same way.

## Editor support

Each example begins with a `yaml-language-server` modeline pointing at the committed schema for its own command in [`packages/schemas`](../packages/schemas), so an editor gives completion, hover documentation, and inline validation with no extra setup. The URL names the command (or `pipeline` for a pipeline file):

```yaml
# yaml-language-server: $schema=https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas/input-config-ancestral.schema.json
```

A JSON config, which cannot carry a comment, uses a top-level `$schema` key with the same URL instead; the loader accepts and ignores it either way. Regenerate the schemas after changing any command's arguments with `treetime schema write --for all -o packages/schemas`.
