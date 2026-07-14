# Output format selection for tree-writing commands

All tree-writing commands (ancestral, clock, timetree, optimize, prune, mugration) produce a fixed set of four files (Newick, Nexus, PhyloGraph JSON, Graphviz DOT) unconditionally. There is no way to select which formats to produce, request additional formats (Auspice, PhyloXML, UShER MAT), choose among Newick annotation styles, or override individual output paths. This proposal introduces a three-tier output format selection system modeled on Nextclade's output architecture.

## Motivation

Three independent problems converge:

1. The shared graph writer needs one format-neutral projection with explicit preservation and loss contracts. Current gaps include [kb/issues/M-core-mutation-representation-and-format-projection-inconsistent.md](../issues/M-core-mutation-representation-and-format-projection-inconsistent.md), [kb/issues/M-io-usher-mat-mutation-loss-is-implicit.md](../issues/M-io-usher-mat-mutation-loss-is-implicit.md), and [kb/issues/N-io-phyloxml-mutation-property-contract-undecided.md](../issues/N-io-phyloxml-mutation-property-contract-undecided.md).
2. Newick annotation styles (plain, BEAST `[&...]`, NHX `[&&NHX:...]`) are a serialization choice that the user cannot control. Plain Newick and annotated Newick are both valid outputs users may need simultaneously.
3. Some commands should produce multiple tree files with different semantics -- timetree should output both time-branch-length and divergence-branch-length trees ([kb/issues/N-io-time-based-branch-lengths-not-implemented.md](../issues/N-io-time-based-branch-lengths-not-implemented.md)), and v0 does this.

A single `--nwk-style` flag cannot address all three. The root problem is that the output system lacks format selection.

## Design

### NWK annotation style as format variant

NWK annotation styles (plain, BEAST-annotated, NHX) are distinct output formats, not a cross-cutting serialization knob. Each variant knows how to serialize itself. Similarly for Nexus, which embeds a Newick string.

This eliminates the `--nwk-style` flag. Instead of "produce Newick, then choose a style," the user selects specific format variants: `nwk`, `nwk-annotated`, `nwk-nhx`.

Rationale: a cross-cutting `--nwk-style` flag applies uniformly to all NWK-family output, preventing "I want the divergence tree as plain NWK but the timetree as annotated NWK." Treating styles as format slots makes every output independently controllable and composes naturally with per-file path overrides.

### Format enum

The `TreeOutputFormat` enum lists all format variants available to analysis commands:

| Variant           | Extension          | Content                                    |
| ----------------- | ------------------ | ------------------------------------------ |
| `nwk`             | `.nwk`             | Plain Newick (no annotations, no comments) |
| `nwk-annotated`   | `.annotated.nwk`   | BEAST-style `[&key=value,...]`             |
| `nwk-nhx`         | `.nhx.nwk`         | NHX `[&&NHX:key=value:...]`                |
| `nexus`           | `.nexus`           | Nexus with plain embedded Newick           |
| `nexus-annotated` | `.annotated.nexus` | Nexus with BEAST-style annotations         |
| `nexus-nhx`       | `.nhx.nexus`       | Nexus with NHX annotations                 |
| `auspice`         | `.auspice.json`    | Auspice v2 JSON                            |
| `phyloxml`        | `.phylo.xml`       | PhyloXML                                   |
| `phyloxml-json`   | `.phylo.json`      | PhyloXML-JSON                              |
| `mat-pb`          | `.mat.pb`          | UShER MAT protobuf                         |
| `mat-json`        | `.mat.json`        | UShER MAT JSON                             |
| `graph-json`      | `.graph.json`      | Internal PhyloGraph JSON                   |
| `dot`             | `.dot`             | Graphviz DOT                               |

Not all variants are available for all commands. The enum is shared; availability is validated per command. See the per-command availability matrix below.

### Per-command format availability

Each cell depends on two factors: whether the format adapter trait is implemented for the command's graph type, and whether the command has the data the format requires. "Works" = trait implemented and data present. "Possible" = data present but trait impl missing. Blank = format not meaningful for that command.

| Format                | ancestral            | optimize             | timetree                | mugration         | clock    | prune    |
| --------------------- | -------------------- | -------------------- | ----------------------- | ----------------- | -------- | -------- |
| nwk/nexus (plain)     | works                | works                | works                   | works             | works    | works    |
| nwk/nexus (annotated) | works (mutations)    | works (mutations)    | works (mutations+dates) | works (states)    |          |          |
| auspice               | possible (mutations) |                      | works                   | possible (traits) |          |          |
| phyloxml              | possible             | possible             | possible                | possible          | possible | possible |
| mat-pb/mat-json       | possible (mutations) | possible (mutations) | possible (mutations)    |                   |          |          |
| graph-json            | works                | works                | works                   | works             | works    | works    |
| dot                   | works                | works                | works                   | works             | works    | works    |

Current annotation data per command (traced from source):

- ancestral: `MutationCommentProvider` emitting `mutations` key (`ancestral/run.rs:155,176`). Empty for parsimony mode
- timetree: `MutationCommentProvider` emitting `mutations` key (`timetree/run.rs:108-109`). `NodeTimetree.nwk_comments()` independently emits `date` key from `self.time` (`payload/timetree.rs:130-135`). Empty when no partitions
- optimize: `MutationCommentProvider` emitting `mutations` key (`optimize/run.rs:62,67`)
- mugration: `DiscreteCommentProvider` emitting `{attribute}` key, e.g. `country` (`mugration/run.rs:75`)
- clock, prune: empty `CommentProviders` -- no annotation data

Format adapter trait implementations that exist today:

- `NodeToNwk`/`EdgeToNwk`, `NodeToGraphviz`/`EdgeToGraphviz`, `Serialize`: all command graph types
- `AuspiceWrite`: timetree only (`TimetreeAuspiceWriter` at `timetree/output/auspice.rs`)
- `PhyloxmlFromGraph`: zero implementations on any analysis command type
- `UsherWrite`: zero implementations on any analysis command type

### Three-tier output control

Follows the Nextclade pattern (`nextclade-cli/src/cli/nextclade_cli.rs`):

Tier 1: bulk directory -- `--output-all` (`-O`): produce all default formats into a directory. Default basename is command-specific (e.g. `timetree`, `annotated_tree`). Replaces the current `--outdir` (`-O`) flag -- the short flag stays the same, the long flag changes to match Nextclade's convention. v1 has not shipped, so no backward compatibility concern.

Tier 2: output selection -- `--output-selection` (`-s`): comma-separated list of output names to include when using `--output-all`. Covers both tree formats and non-tree outputs in a single enum (matching Nextclade's `--output-selection` pattern). Default tree formats: `nwk,nexus`. This intentionally drops `graph-json` and `dot` from the current default set (`write_graph_files_with_options` at `treetime-io/src/graph.rs:18-53` unconditionally writes all four). Example: `--output-selection=nwk,nwk-annotated,auspice,augur-node-data`.

Tier 3: per-file path overrides -- one flag per format variant: `--output-tree-nwk=<path>`, `--output-tree-nwk-annotated=<path>`, `--output-tree-auspice=<path>`, etc. Each takes an `Option<PathBuf>`. A per-file flag produces that format regardless of Tier 2 selection, at the specified path.

Resolution logic (same as Nextclade's `get_or_insert` pattern):

1. If `--output-all` is set, compute default paths: `<dir>/<basename>.<extension>` for each format in Tier 2 selection
2. Per-file flags override or supplement: if the flag is already `Some`, default is not applied (`get_or_insert` semantics)
3. Per-file flags outside the Tier 2 selection still produce their format (explicit request always honored)
4. No output flags at all -> error (at least one output required)

### Command-specific tree identities

Commands that produce multiple tree outputs (different branch lengths, different topologies) declare named tree identities. Each identity gets the full format treatment.

| Command   | Tree identities               | Branch length semantics        |
| --------- | ----------------------------- | ------------------------------ |
| timetree  | `timetree`, `divergence-tree` | Time units, substitutions/site |
| ancestral | `annotated-tree`              | Substitutions/site             |
| optimize  | `annotated-tree`              | Substitutions/site             |
| mugration | `annotated-tree`              | Substitutions/site             |
| clock     | `rerooted-tree`               | Substitutions/site             |
| prune     | `pruned-tree`                 | Substitutions/site             |

Per-file flags for multi-tree commands use the tree identity as prefix: `--output-timetree-nwk`, `--output-divergence-tree-nwk-annotated`. Single-tree commands use the generic `--output-tree-*` prefix.

### Non-tree outputs

Flat per-file flags, no format axis. `--output-all` produces all non-tree outputs that the command supports alongside the selected tree formats. Per-file flags override or suppress individual non-tree outputs.

- `--output-augur-node-data=<path>` (ancestral, timetree)
- `--output-fasta=<path>` (ancestral)
- `--output-clock-model=<path>` (timetree, clock)
- `--output-traits-csv=<path>` (mugration)
- `--output-gtr=<path>` (ancestral, optimize, timetree, mugration)
- `--output-confidence=<path>` (timetree, mugration)

### Compression

Transparent compression from output path extension (`.gz`, `.bz2`, `.xz`, `.zst`), inherited from the existing I/O layer. `-` writes to stdout (uncompressed).

## Writer dispatch

The current NWK writer hardcodes BEAST-style annotations at `treetime-io/src/nwk.rs:255`:

```rust
.map(|(key, val)| format!("[&{key}=\"{val}\"]"))
```

`CommentProviders` produces dialect-agnostic `BTreeMap<String, String>` data per node. The format variant determines how that data is serialized: BEAST `[&key="value"]`, NHX `[&&NHX:key=value:...]`, or suppressed (plain). This single format string becomes a dispatch on the output format variant.

Annotation data reaches the writer from two independent sources that merge at `nwk.rs:239`:

- `NodeToNwk::nwk_comments()` on the node type (e.g. `NodeTimetree` emits `date` at `payload/timetree.rs:130-135`)
- `CommentProviders` injected by the command (e.g. `MutationCommentProvider` emits `mutations`)

Both produce `BTreeMap<String, String>`. The format variant dispatch applies after merging, at the serialization step.

## Impact

- Resolves [kb/issues/N-io-time-based-branch-lengths-not-implemented.md](../issues/N-io-time-based-branch-lengths-not-implemented.md) -- timetree outputs both tree identities
- Intentional default change: `graph-json` and `dot` no longer produced by default (currently produced unconditionally by `write_graph_files_with_options`)

## Validation

- Each `TreeOutputFormat` variant round-trips every value in its declared vocabulary. Mutation round trips remain conditional on resolving the shared mutation vocabulary and each writer's explicit unsupported-state policy; schema-level serialization alone does not establish semantic round-trip correctness.
- Tier 1/2/3 resolution logic: unit tests for flag precedence, format selection filtering, default path generation
- Per-command: verify available format set matches command's adapter traits
- v0 parity: timetree produces both `timetree.nwk` and `divergence_tree.nwk` (matching v0's `export_sequences_and_tree`)

## Related

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) -- format adapter architecture
- [kb/reports/newick-annotation-dialects.md](../reports/newick-annotation-dialects.md) -- NWK dialect grammars and tool interop
- [kb/proposals/unified-input-format-support.md](unified-input-format-support.md) -- input-side counterpart (analysis commands accept any format)
- [kb/tickets/io-timetree-divergence-tree-output.md](../tickets/io-timetree-divergence-tree-output.md) -- timetree divergence-tree output
