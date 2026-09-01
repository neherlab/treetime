# Example configs and pipelines

This directory holds the example datasets and, next to each one, ready-to-run configuration files for TreeTime. Every command reads its full argument object from a file with `--config`, and the `pipeline` command runs an ordered chain of commands from a single file. The examples double as living documentation: each file is a saved, reviewed, version-controlled invocation you can copy and adapt.

All examples are run from the repository root and write their outputs under `tmp/` (git-ignored). Input paths inside the files are relative to the repository root.

## Per-command configs

A per-command config is exactly the command's argument object. Its keys map one-to-one to the CLI flags: `--alignment a.fasta` is `alignment: [a.fasta]`, `--clock-rate 0.003` is `clock_rate: 0.003`, `--output-all dir` is `output_all: dir`. An explicit flag on the command line still overrides the file (flag > config > default).

| File                                                                     | Command     | Demonstrates                                                |
| ------------------------------------------------------------------------ | ----------- | ----------------------------------------------------------- |
| [`ebola/20/ancestral-parsimony.yaml`](ebola/20/ancestral-parsimony.yaml) | `ancestral` | Fitch parsimony: fast, model-free ancestral states          |
| [`zika/20/ancestral.yaml`](zika/20/ancestral.yaml)                       | `ancestral` | Marginal reconstruction with an inferred GTR model          |
| [`sc2/2844/ancestral.yaml`](sc2/2844/ancestral.yaml)                     | `ancestral` | Marginal reconstruction on 2844 SARS-CoV-2 genomes          |
| [`dengue/100/clock.yaml`](dengue/100/clock.yaml)                         | `clock`     | Root-to-tip regression with a non-default id column         |
| [`flu/h3n2/20/optimize.yaml`](flu/h3n2/20/optimize.yaml)                 | `optimize`  | Branch-length optimization and min-divergence rerooting     |
| [`mpox/clade-ii/1000/optimize.yaml`](mpox/clade-ii/1000/optimize.yaml)   | `optimize`  | The mpox `augur refine` replacement (mutations, no indels)  |
| [`rsv/a/100/prune.yaml`](rsv/a/100/prune.yaml)                           | `prune`     | Dropping short, empty, and shared-mutation branches         |
| [`zika/20/mugration.yaml`](zika/20/mugration.yaml)                       | `mugration` | Discrete ancestral geography over sampling country          |
| [`flu/h3n2/200/timetree.yaml`](flu/h3n2/200/timetree.yaml)               | `timetree`  | Dated tree with a fixed clock rate and confidence intervals |

Run one with, for example:

```bash
treetime ancestral --config data/zika/20/ancestral.yaml
```

## Pipelines

A pipeline file lists named steps, each an analysis command with its arguments. Steps run sequentially in one process. A step reads an earlier step's output with `{{ steps.<name>.outputs.<selection> }}`, so intermediate paths are never hand-written. `vars` names shared values once; a top-level `output_all` gives each step its own `<dir>/<step>/` output directory.

| File                                                                   | Chain                                   | Demonstrates                                                     |
| ---------------------------------------------------------------------- | --------------------------------------- | ---------------------------------------------------------------- |
| [`zika/20/pipeline.yaml`](zika/20/pipeline.yaml)                       | `timetree` -> `mugration`               | Phylodynamics and phylogeography end to end                      |
| [`flu/h3n2/200/pipeline.yaml`](flu/h3n2/200/pipeline.yaml)             | `optimize` -> `timetree`                | A shared, typed clock rate defined once in `vars`                |
| [`flu/h3n2/500/pipeline.yaml`](flu/h3n2/500/pipeline.yaml)             | `optimize` -> `ancestral` -> `timetree` | A seasonal-influenza build with a fixed clock (500 sequences)    |
| [`ebola/362/pipeline.yaml`](ebola/362/pipeline.yaml)                   | `optimize` -> `ancestral` -> `timetree` | A three-step build on the 2014 Ebola epidemic dataset            |
| [`sc2/2844/pipeline.yaml`](sc2/2844/pipeline.yaml)                     | `optimize` -> `ancestral`               | A date-free build on 2844 SARS-CoV-2 genomes                     |
| [`mpox/clade-ii/1000/pipeline.yaml`](mpox/clade-ii/1000/pipeline.yaml) | `optimize` -> `ancestral` -> `timetree` | The full Nextstrain mpox TreeTime chain (fixed clock, mutations) |

Run a pipeline, preview its plan, or run only part of it:

```bash
treetime pipeline --config data/zika/20/pipeline.yaml
treetime pipeline --config data/zika/20/pipeline.yaml --check
treetime pipeline --config data/zika/20/pipeline.yaml --steps mugration
```

`--check` resolves the whole plan (each step's inputs, its concrete output files, and every chained path) and prints it without running anything. `--steps` runs a named subset in list order; a referenced upstream step that is not selected must already have its outputs on disk. When a step fails, the error names the completed steps and prints the exact `--steps` value that resumes the rest.

## Why a config or pipeline instead of a shell command

- One file, not one long flag string. The invocation is reviewable in a pull request, diffable across runs, and stored next to the data it analyzes.
- No repetition. `vars` holds the dataset directory or the clock rate once; every step references it. A whole-value reference such as `{{ vars.clock_rate }}` keeps its JSON type, so a number stays a number.
- No path bookkeeping. `{{ steps.optimize.outputs.nwk }}` resolves to the optimized tree the first step actually writes. Rename an output directory in one place and the chain follows.
- Errors before computation. `--check` resolves and type-checks the whole plan up front. A misspelled key is rejected with a suggestion (`unknown field opt_methd`; `unknown top-level key step; did you mean steps?`) rather than running silently or halfway.
- One process. Steps share the thread pool; there is no shell loop, no temporary files to name, and no second tool to install.

## Why a pipeline instead of Snakemake for the TreeTime part

- No workflow engine to stand up. A Snakemake build needs a `Snakefile`, a `config.yaml`, wildcards, and one rule per tool. The pipeline file needs only the TreeTime binary.
- Chaining by output identity, not by string. A Snakemake rule wires steps with literal `input:`/`output:` path strings that must match each tool's real filenames. The pipeline references an output by its selection tag (`nwk`, `auspice`, ...) and resolves the concrete path itself.
- Validation is static and typed. `--check` and the JSON-schema-backed argument objects catch a bad parameter or a dangling step reference before anything runs; a Snakemake rule only surfaces a bad TreeTime parameter when that rule executes.
- Reproducible and portable. The config is the exact serialized argument shape, with one dependency, committed beside the data.

This is not a replacement for a full multi-sample workflow manager: incremental rebuilds, cluster scheduling, and non-TreeTime rules (subsampling, alignment, tree building, export) still belong to Snakemake and Augur. The pipeline replaces the TreeTime glue inside such a build. On the Nextstrain mpox `WIP/new_tt` phylogenetic workflow, for instance, three Snakemake rules drive TreeTime: `treetime optimize` replaces `augur refine` (branch-length optimization), `treetime ancestral` folds in `augur translate`, and `treetime timetree` dates the tree. The [`mpox/clade-ii/1000/pipeline.yaml`](mpox/clade-ii/1000/pipeline.yaml) example runs that chain as three steps in one file.

## Editor support

Both file kinds are plain JSON or YAML (one parser reads both). Generate the JSON schemas that back completion and validation with:

```bash
treetime schema --for all -o tmp/schemas
```

A pipeline file may carry a top-level `$schema` key pointing at `pipeline.schema.json`; the loader ignores it and editors use it. A per-command config has no `$schema` key (the strict commands reject unknown keys), so map its schema in the editor by filename instead, for example with the YAML extension's `yaml.schemas` setting keyed on `optimize.schema.json`, `timetree.schema.json`, and so on.
