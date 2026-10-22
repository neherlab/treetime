# CLI help text defects, inconsistencies, and UX violations

> [!IMPORTANT]
> **Decision required.** Each behavioral item needs its own parse and runtime contract before implementation, and the split of this inventory is unapproved: focused issues per behavior contract (O1), or one repository-wide CLI cleanup (O2). See "Potential solutions".

Systematic audit of `--help` output across all commands reveals defects affecting human scientists, automation clients, and workflow engines such as Snakemake and Nextflow.

## Defects

### D7: `--alphabet` help describes detection that does not exist

The help of `--alphabet` says that the alphabet is detected from the sequence content when the flag is omitted, and the doc comment of the argument struct names `detect_alphabet` [packages/app-commands/src/commands/shared/alphabet.rs#L8-L11](../../packages/app-commands/src/commands/shared/alphabet.rs#L8-L11) [packages/app-commands/src/commands/shared/alphabet.rs#L22-L25](../../packages/app-commands/src/commands/shared/alphabet.rs#L22-L25). No such function exists: every caller uses `alphabet_name().unwrap_or_default()`, which gives `nuc`. The same text reaches the generated CLI reference, the JSON schemas, and the TypeScript types.

## Inconsistencies

### I6: Capitalization inconsistent across all commands

Some descriptions of `timetree` start lowercase ("use an autocorrelated"), others uppercase ("If set to 'input'", "Method used for", "Length of the sequence").

## UX / ergonomics

### U2: `--optimizer` renders LaTeX in terminal

Command: `optimize`. The description and possible value descriptions contain `$\sqrt{t}$`, `$s = \sqrt{t}$`, `$O(1/t^2)$` which display as raw LaTeX in terminal output. Should be `sqrt(t)`, `s = sqrt(t)`, `O(1/t^2)`.

### U4: `--relax` uses unusual multi-value syntax

Command: `timetree`. Takes two positional-style values (`--relax 1.0 0.5`). Most CLIs use separate flags or comma-separated values.

## Style

### S3: `mugration --weights` grammar error

"probabilities of that a randomly sampled sequence..." should be "probability that a randomly sampled sequence...".

### S4: `--alphabet aa-no-stop` unexplained

No description of what `aa-no-stop` means or when to choose it over `aa`.

### S5: Various descriptions reference v0 internals

`--smooth-initial-pi`, `--filter-uninformative-root` in `mugration` reference "TreeTime v0" behavior and `infer_gtr regularization`.

### S6: `--sampling-bias-correction` vague

No indication of expected range, units, or typical values.

## Workflow/automation

### W1: Output flags appear optional but at least one is required at runtime

Usage lines show all output flags as optional (`[OPTIONS]`), but the output plan returns `MissingOutputs::NoOutputs` when no output destination is provided [packages/app-output/src/output_plan.rs#L91-L93](../../packages/app-output/src/output_plan.rs#L91-L93). Users following the usage pattern get an argument error. Workflow generators cannot infer that `-O`/`--output-all` or a per-file flag is required.

### W2: `--output-selection` accepts values that the command rejects at runtime

The `command_outputs!` macro adds the full tree-format set (`nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`) to the selection enum of every command [packages/app-commands/src/commands/shared/output_args.rs#L22-L38](../../packages/app-commands/src/commands/shared/output_args.rs#L22-L38), while a separate availability matrix rejects unsupported values later. Per-file flags have the same drift when their writer is absent from the command matrix.

Generate each command-specific parseable enum and visible per-file flags from its actual format set, then add parse-level rejection tests.

### W6: No usage examples in any scientific command

Core commands (`timetree`, `ancestral`, `clock`, `mugration`, `optimize`, `prune`) have no canonical command-line examples showing required input combinations. The help texts with examples are `completions` and the `--model-params` flag (which shows `--model hky85 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`).

### W7: `mugration` metadata format examples are unreadable inline prose

`--metadata` description: "CSV or TSV file with discrete characters. #name,country,continent taxon1,micronesia,oceania ..." -- column structure is not visually separated from data rows.

### W8: `schema` and `help-markdown` mixed with scientific commands

Root help lists tooling commands at the same level as analysis commands with no grouping or labeling.

## Grammar

### G2: "gaussian" should be "Gaussian"

[packages/app-commands/src/commands/timetree/args.rs#L215](../../packages/app-commands/src/commands/timetree/args.rs#L215)

### G4: Root `ancestral` summary is a full paragraph; other summaries are one-liners

[packages/app-cli/src/cli/treetime_cli.rs#L126](../../packages/app-cli/src/cli/treetime_cli.rs#L126): uneven length makes root help hard to scan. The `clock` and `homoplasy` summaries are also several sentences long.

## Affected code

- CLI arg definitions: `packages/app-commands/src/commands/*/args.rs`
- Shared arg structs: `packages/app-commands/src/commands/shared/*.rs`
- Top-level command descriptions: `packages/app-cli/src/cli/treetime_cli.rs`
- Method enum: `packages/treetime/src/ancestral/params.rs`

## Potential solutions

- O1. Extract one issue per independent parser or behavior contract, then group only purely textual corrections that share the same help generator and validation snapshot.
- O2. Keep one repository-wide CLI cleanup. This couples parser behavior, dormant feature decisions, generated reference output, and prose edits into an implementation that cannot be validated as one contract.

## Recommendation

Use O1. Parser behavior such as dense/sparse mode, required outputs, and per-command format availability needs focused issues. Pure spelling, capitalization, grammar, and layout corrections may share one documentation change after behavior-dependent items are removed. This inventory is not executable as one change.

## Readiness

Each behavioral item must first select its parse/runtime contract, and already-owned unused or unimplemented flags remain with their domain issues.

## Related issues

- [M-cli-flags-parsed-but-ignored.md](M-cli-flags-parsed-but-ignored.md): unused flags
