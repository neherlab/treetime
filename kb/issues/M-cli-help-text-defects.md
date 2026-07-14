# CLI help text defects, inconsistencies, and UX violations

Systematic audit of `--help` output across all commands reveals defects affecting human scientists, automation clients, and workflow engines such as Snakemake and Nextflow.

## Defects

### D1: `--method-anc=joint` accepted by parser, errors at runtime

Commands: `timetree`, `ancestral`, `clock`, `homoplasy`.

`[possible values: marginal, parsimony, joint]` lists `joint`, but the pipeline rejects it: "Joint ancestral reconstruction has been removed. Available methods: marginal, parsimony" (`ancestral/pipeline.rs:231`). Users and automation discover this only at runtime, after setup completes.

Fix: remove `joint` from the enum or hide it with a clap `hide` attribute and add a deprecation note.

### D2: `--tree` description promises fallback on required-tree commands

Commands: `optimize`, `ancestral`.

Description says "If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml" but `--tree` is **required** (`PathBuf`, shown as `--tree <TREE>` in usage line). The fallback text is unreachable.

Fix: use a different description for commands where `--tree` is required. `prune` already does this correctly ("Name of file containing the tree in newick, nexus, or phylip format" with no fallback clause).

Related: [H-timetree-tree-inference-unimplemented.md](H-timetree-tree-inference-unimplemented.md) covers the `timetree`/`clock` case where `--tree` is optional but tree inference is not implemented.

### D3: `--prune-short` in `clock` has empty description and type mismatch

`clock` command: `--prune-short` renders with no description at all (no doc comment on `commands/clock/args.rs:76`). It is a `bool` flag. In `prune`, the same name takes an `Option<f64>` threshold (`--prune-short <THRESHOLD>`). Same flag name, different types, different semantics.

Related: [M-clock-dead-cli-arguments.md](M-clock-dead-cli-arguments.md) tracks that `--prune-short` in `clock` is parsed but never wired.

### D4: `clock` command description has `--keep_root` typo

Top-level `clock` description: "unless run with --keep_root" (underscore). The actual flag is `--keep-root` (hyphen). Copy-pasting from the description will fail at parse time.

### D5: `--model-params` references Python source files

Commands: `timetree`, `optimize`, `ancestral`, `clock`, `homoplasy`.

Description: "See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py". These are v0 Python paths that do not exist in the Rust codebase. Users of the Rust binary cannot find these files.

Fix: reference the v1 Rust module paths or link to documentation.

### D7: `ancestral` command description hardcodes output filenames

Description: "The output consists of a file 'ancestral.fasta' with ancestral sequences and a tree 'ancestral.nexus'". Actual output is controlled by `--output-all` and per-file flags. These specific filenames are not guaranteed. (Partially fixed: old `annotated_tree.nexus` reference updated to `ancestral.nexus` to match S4 stem change.)

### D8: `--n-iqd` in `timetree` is unused

Defined at `commands/timetree/args.rs:184` but never read or used anywhere. The `--clock-filter` flag already controls interquartile-based outlier detection. The relationship between the two is unexplained.

Related: [N-timetree-dead-cli-flags.md](N-timetree-dead-cli-flags.md) lists `--n-iqd`.

### D9: `homoplasy` and `arg` show full help for unimplemented commands

`homoplasy` displays a complete options page but `run_homoplasy` returns "not yet implemented in v1". `arg` shows global flags only with no indication of being unimplemented.

Related: [H-homoplasy-command-unimplemented.md](H-homoplasy-command-unimplemented.md).

## Inconsistencies

### I1: `--metadata` short flag and alias differ across commands

`timetree`/`clock`: `-d --metadata` with alias `--dates`. `mugration`: `-s --metadata` with alias `--states`. Same long flag name, different short flags.

### I2: `--clock-filter` description differs between `timetree` and `clock`

"inter-quartile" (hyphenated) in `timetree` vs "interquartile" (one word) in `clock`. Both include redundant prose "Default=3.0" alongside the `[default: 3.0]` metadata.

### I3: `--dense` description differs across four commands

`timetree`: "store full probability distributions". `optimize`: "stores full probability vectors ... Dense is more accurate..." (extra paragraph). `ancestral`/`homoplasy`: "stores full probability vectors at each position" (extra GTR paragraph).

### I4: `--vcf-reference` description has two phrasings

`timetree`/`clock`: "Only for vcf input: fasta file of the sequence the VCF was mapped to". `ancestral`/`homoplasy`: "FASTA file of the sequence the VCF was mapped to (only for vcf input)".

### I5: `--max-iter` default and description differ

`timetree`: default 2, "maximal number" (lowercase). `optimize`: default 10, "Maximum number" (capitalized). Different defaults are intentional but the descriptions diverge without reason.

### ~~I6: `--confidence` has incompatible semantics across commands~~ (resolved)

Resolved by splitting into `--output-confidence-tsv` (timetree) and `--output-confidence-csv` (mugration, alias `--confidence`).

### I7: Capitalization inconsistent across all commands

Some descriptions start lowercase ("don't reroot the tree", "ignore tips", "excess variance", "maximal number", "use an autocorrelated", "rescale branch lengths"), others uppercase ("If set to 'input'", "Method used for", "Length of the sequence").

### ~~I8: `--output-selection` shows full superset regardless of command~~ (resolved)

Resolved by per-command selection enums (`AncestralOutputSelection`, `TimetreeOutputSelection`, etc.) that restrict `--output-selection` to only the outputs the command produces.

### I9: `--branch-length-mode` description differs subtly

`timetree`: "Branch lengths optimized by treetime are only accurate at short evolutionary distances". `clock`: "Note that branch lengths optimized by treetime are only accurate at short evolutionary distances".

## UX / ergonomics

### U1: `--dense` takes `true`/`false` string instead of boolean flag

Commands: `timetree`, `optimize`, `ancestral`, `homoplasy`. `[possible values: true, false]` requires `--dense=true` or `--dense true`. Idiomatic CLI would be `--dense`/`--no-dense` flag pair or a tri-state with auto-detection.

### U2: `--opt-method` renders LaTeX in terminal

Command: `optimize`. Possible value descriptions contain `$t$`, `$\sqrt{t}$`, `$\ln(t)$` which display as raw LaTeX in terminal output. Should be `t`, `sqrt(t)`, `ln(t)`.

### U4: `--relax` uses unusual multi-value syntax

Command: `timetree`. Takes two positional-style values (`--relax 1.0 0.5`). Most CLIs use separate flags or comma-separated values.

### U5: `--detailed` in `homoplasy` takes opaque string value

`--detailed <DETAILED>` with description "generate a more detailed report". No documentation of valid values.

### U6: `--pc` in `mugration` default not shown in metadata

Description says "Default: 1.0" but no `[default: 1.0]` metadata (field is `Option<f64>`). Automation parsers miss it.

### U7: `--keep-polytomies`/`--resolve-polytomies` no conflict guard

Command: `timetree`. Both are `bool` flags with no `conflicts_with`. Both can be set simultaneously.

Related: [N-timetree-polytomy-flags-no-conflict.md](N-timetree-polytomy-flags-no-conflict.md).

### U8: Global flags on non-analysis commands

`--jobs`, `--verbosity`, `--no-progress` appear on `completions`, `schema`, `help-markdown`, `arg` where they serve no purpose.

## Style

### S1: Format name capitalization

"newick, nexus, or phylip" should be "Newick, Nexus, or PHYLIP" consistently.

### S2: "compressed fasta" should be "compressed FASTA"

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

Usage lines show all output flags as optional (`[OPTIONS]`), but `output.rs:562` errors when no output destination is provided. Users following the usage pattern get a runtime error after input loading completes. Workflow generators cannot infer that `-O`/`--output-all` or a per-file flag is required.

### W2: `--output-selection` accepts values that the command rejects at runtime

The macro still adds the full tree-format superset to every command, while a separate availability matrix rejects unsupported values later. For example, `prune --help` advertises Auspice and MAT selections that prune cannot produce. Per-file flags have the same drift when their writer is absent from the command matrix.

Generate each command’s parseable enum and visible per-file flags from its actual format set, then add parse-level rejection tests.

### W3: `--alignment` says "multiple FASTA files" but does not show repeat syntax

Help renders `-a, --alignment <FILEPATH>` with singular value name. No indication whether multiple files use repeated flags (`--alignment a --alignment b`), space-separated values, or comma-separated values.

### W4: Stdin fallback text appears on commands where it is irrelevant

Shared alignment help says "If no input files provided, the plain fasta input is read from standard input (stdin)." This appears on `clock` (metadata/tree-based, not sequence-based) and `prune` (sequence only needed for `--prune-empty`).

### W5: Root help links to Python v0 documentation

Root `--help` links to `https://treetime.readthedocs.io/en/stable/` and the Python TreeTime publication. Readers can incorrectly infer v1 behavior from documentation for the v0 binary.

### W6: No usage examples in any scientific command

Only `completions` includes an example. Core commands (`timetree`, `ancestral`, `clock`, `mugration`, `optimize`, `prune`) have no canonical command-line examples showing required input combinations.

### W7: `mugration` metadata format examples are unreadable inline prose

`--metadata` description: "CSV or TSV file with discrete characters. #name,country,continent taxon1,micronesia,oceania ..." -- column structure is not visually separated from data rows.

### W8: `schema` and `help-markdown` mixed with scientific commands

Root help lists tooling commands at the same level as analysis commands with no grouping or labeling.

## Grammar

### G1: "If there's multiple input files" should be "If there are multiple input files"

`commands/shared/alignment.rs:21`

### G2: "gaussian" should be "Gaussian"

`commands/timetree/args.rs:102`

### G3: "Higher numbers results" should be "Higher numbers result"

`commands/mugration/args.rs:44`

### G4: Root `ancestral` summary is a full paragraph; other summaries are one-liners

`app-cli/src/cli/treetime_cli.rs:73` -- uneven length makes root help hard to scan.

## Affected code

- CLI arg definitions: `packages/treetime/src/commands/*/args.rs`
- Shared arg structs: `packages/treetime/src/commands/shared/*.rs`
- Top-level command descriptions: `packages/app-cli/src/cli/treetime_cli.rs`
- Method enum: `packages/treetime/src/ancestral/params.rs`

## Potential solutions

- O1. Extract one issue per independent parser or behavior contract, then group only purely textual corrections that share the same help generator and validation snapshot.
- O2. Keep one repository-wide CLI cleanup. This couples parser behavior, dormant feature decisions, generated reference output, and prose edits into an implementation that cannot be validated as one contract.

## Recommendation

Use O1. Parser behavior such as `joint`, dense/sparse mode, required outputs, and per-command format availability needs focused issues and tickets. Pure spelling, capitalization, grammar, and layout corrections may share one documentation ticket after behavior-dependent items are removed. This inventory has no executable omnibus ticket.

## Ticket readiness

No aggregate ticket is ready. Each behavioral item must first select its parse/runtime contract, and already-owned dead or unimplemented flags remain with their domain issues.

## Related issues

- [N-timetree-dead-cli-flags.md](N-timetree-dead-cli-flags.md): unused flags in timetree (overlaps D8)
- [M-clock-dead-cli-arguments.md](M-clock-dead-cli-arguments.md): unused flags in clock (overlaps D3)
- [H-timetree-tree-inference-unimplemented.md](H-timetree-tree-inference-unimplemented.md): tree inference fallback not implemented (overlaps D2)
- [H-homoplasy-command-unimplemented.md](H-homoplasy-command-unimplemented.md): homoplasy unimplemented (overlaps D9)
- [N-timetree-polytomy-flags-no-conflict.md](N-timetree-polytomy-flags-no-conflict.md): polytomy flag conflict (overlaps U7)
- [M-timetree-method-anc-ignored.md](M-timetree-method-anc-ignored.md): `--method-anc` dead in timetree (related to D1)
