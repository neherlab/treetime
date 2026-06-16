# Newick branch length precision

The `treetime-io` Newick writer defaults to 3 significant digits for branch lengths, the lowest precision of any phylogenetic tool surveyed. This silently destroys information on every write cycle. Branch lengths of `0.123456789` become `0.123` -- a 0.4% relative error that compounds across downstream analyses reading v1 output. v0 uses 7 decimal places for divergence trees; the Nextstrain ecosystem (Augur) uses 8 decimal places; four major tools default to full f64 precision.

## Problem

`fn format_weight()` [[src](https://github.com/neherlab/treetime/blob/d9e5936b/packages/treetime-io/src/nwk.rs#L295-L304)] applies a hardcoded `Some(3)` fallback when no precision is configured:

```rust
pub fn format_weight(weight: f64, options: &NwkWriteOptions) -> String {
  float_to_digits(
    weight,
    options.weight_significant_digits.or(Some(3)),  // <- hardcoded 3
    options.weight_decimal_digits,
  )
}
```

`struct NwkWriteOptions` [[src](https://github.com/neherlab/treetime/blob/d9e5936b/packages/treetime-io/src/nwk.rs#L114-L123)] defaults both `weight_significant_digits` and `weight_decimal_digits` to `None`. All command output flows through `fn write_tree_outputs()` [[src](https://github.com/neherlab/treetime/blob/f5710471/packages/treetime-io/src/graph.rs#L57-L62)] which constructs `NwkWriteOptions::default()`, leaving both fields as `None` and triggering the 3-sigfig fallback. The Graphviz DOT writer [[src](https://github.com/neherlab/treetime/blob/b6c643bd/packages/treetime-io/src/graphviz.rs#L189)] also calls `format_weight` with `NwkWriteOptions::default()`, so it has the same truncation.

A second instance of the same default exists in `fn float_to_digits()` [[src](https://github.com/neherlab/treetime/blob/0535f431/packages/treetime-utils/src/fmt/float.rs#L73-L94)] at line 82: when neither constraint is specified, it also falls back to 3 significant digits. This affects any caller of the general-purpose float formatter.

The standalone `util-newick` crate's `fn format_float()` [[src](https://github.com/neherlab/treetime/blob/d9e5936b/packages/util-newick/src/write.rs)] does not have this problem -- it passes `None`/`None` through to `pretty_dtoa` which then uses full precision.

## Ecosystem survey

All tools inspected from source at commit SHAs listed below. Claims in the original issue file ([kb/issues/N-io-nwk-writer-3-sigfig-default-truncates-precision.md](../issues/N-io-nwk-writer-3-sigfig-default-truncates-precision.md)) that DendroPy uses `%.10e` by default are incorrect -- DendroPy defaults to Python's `str(float)` (full precision).

### Full-precision defaults

| Tool     | Mechanism                                                                                   | Configurable                            | Source                                                                                                                                                                                                                    |
| -------- | ------------------------------------------------------------------------------------------- | --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| gotree   | `strconv.FormatFloat(length, 'f', -1, 64)` -- shortest exact round-trip representation      | no                                      | [[src](https://github.com/evolbioinfo/gotree/blob/83cca67f/tree/node.go#L259)]                                                                                                                                            |
| BEAST2   | `StringBuilder.append(double)` -> Java `Double.toString()` -- shortest exact representation | yes (`dp` parameter, default -1 = full) | [[src](https://github.com/CompEvol/beast2/blob/9321c88d/src/beast/base/evolution/tree/Node.java#L404)] [[src](https://github.com/CompEvol/beast2/blob/9321c88d/src/beast/base/evolution/TreeWithMetaDataLogger.java#L29)] |
| BEAST1   | `PrintStream.print(double)` -> `Double.toString()`                                          | yes (`setNumberFormat()`)               | [[src](https://github.com/beast-dev/beast-mcmc/blob/24826333/src/dr/app/tools/NexusExporter.java#L296-L302)]                                                                                                              |
| DendroPy | `str(float)` -- Python repr, full precision                                                 | yes (`real_value_format_specifier`)     | [[src](https://github.com/jeetsukumaran/DendroPy/blob/ffffe49f/src/dendropy/dataio/newickwriter.py#L189)]                                                                                                                 |

### Fixed-precision defaults

| Tool             | Default                   | Type                                                                                | Configurable                                       | Source                                                                                                                                                                        |
| ---------------- | ------------------------- | ----------------------------------------------------------------------------------- | -------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| IQ-TREE          | 10 decimal places         | `std::fixed` with `precision(10)`                                                   | adaptive (scales with `min_branch_length`)         | [[src](https://github.com/iqtree/iqtree2/blob/6776a95f/tree/mtree.cpp#L430-L447)]                                                                                             |
| Augur            | 8 decimal places          | `%1.8f`                                                                             | no                                                 | [[src](https://github.com/nextstrain/augur/blob/024292af/augur/tree.py#L564)] [[src](https://github.com/nextstrain/augur/blob/024292af/augur/refine.py#L438)]                 |
| TreeTime v0      | 7 decimal / 9 significant | `%1.7f` (short aln) or `%1.9e` (long aln); timetree falls back to Biopython `%1.5f` | no                                                 | [src](../../packages/legacy/treetime/treetime/CLI_io.py#L209)                                                                    |
| RAxML-NG         | 6 decimal places          | `%.*lf` with `RAXML_DEFAULT_PRECISION = 6`                                          | yes (`--precision N`); adaptive with `--brlen-min` | [[src](https://github.com/amkozlov/raxml-ng/blob/d850b2cf/src/io/NewickStream.cpp#L10-L21)] [[src](https://github.com/amkozlov/raxml-ng/blob/d850b2cf/src/constants.hpp#L44)] |
| MrBayes          | 6 decimal (scientific)    | `%.*le`                                                                             | yes (`set precision=N`, range 3-15)                | [[src](https://github.com/NBISweden/MrBayes/blob/bb09fffb/src/utils.c#L1041-L1051)]                                                                                           |
| ETE4             | 6 significant digits      | `%g`                                                                                | yes (parser dict)                                  | [[src](https://github.com/etetoolkit/ete/blob/9562dfb6/ete4/parser/newick.pyx#L111)]                                                                                          |
| Biopython        | 5 decimal places          | `%1.5f`                                                                             | yes (`format_branch_length` kwarg)                 | [[src](https://github.com/biopython/biopython/blob/75bc2aea/Bio/Phylo/NewickIO.py#L281)]                                                                                      |
| **v1 (current)** | **3 significant digits**  | `pretty_dtoa` with `max_significant_digits(3)`                                      | fields exist but not wired                         | [[src](https://github.com/neherlab/treetime/blob/d9e5936b/packages/treetime-io/src/nwk.rs#L295-L304)]                                                                         |

v1 is the least-precise writer in the survey by a substantial margin.

### Precision comparison on representative branch lengths

| Value       | v1 (3 sig) | v0 (7 dec) | Augur (8 dec) | Full        |
| ----------- | ---------- | ---------- | ------------- | ----------- |
| 0.123456789 | 0.123      | 0.1234568  | 0.12345679    | 0.123456789 |
| 0.00123456  | 0.00123    | 0.0012346  | 0.00123456    | 0.00123456  |
| 0.0001234   | 0.000123   | 0.0001234  | 0.00012340    | 0.0001234   |
| 0.99876543  | 0.999      | 0.9987654  | 0.99876543    | 0.99876543  |

## Design decisions

Three independent axes. Each can be decided separately.

### Axis 1: default precision level

What precision does `treetime-io` use when no explicit constraint is set?

#### A. Full precision (shortest round-trip, Ryu-based)

`pretty_dtoa` with no digit limit produces the shortest decimal string that round-trips to the same f64 value (typically 15-17 significant digits, trailing zeros trimmed). This is what `util-newick`'s writer already does.

- Pro: lossless -- write-read-write cycles preserve values exactly
- Pro: simplest implementation -- remove the fallback, let `pretty_dtoa` do its default
- Pro: aligns with gotree, BEAST1/2, DendroPy (4 of 10 surveyed tools)
- Pro: aligns with `util-newick` writer already in this codebase, eliminating inconsistency between the two Newick writers
- Pro: overshoots v0 rather than undershooting -- strictly more information preserved
- Con: output files slightly larger (e.g. `0.00123456789012345` vs `0.00123`)
- Con: output differs visually from v0 (more digits than `%1.7f`)

#### B. Match v0 precision (7 decimal / 9 significant, adaptive)

Replicate v0's `%1.7f` (short alignments) / `%1.9e` (long alignments) logic.

- Pro: golden master tests pass without tolerance adjustments
- Pro: output visually identical to v0
- Con: requires alignment-length-dependent formatting logic in the I/O layer, coupling serialization to analysis state
- Con: still lossy for branch lengths near 1e-7 (7 decimal places truncates `1.23e-8` to `0.0000000`)
- Con: replicates a Biopython-specific pattern (`%1.7f` / `%1.9e`) with no scientific justification -- the digit count was chosen for Biopython's `format_branch_length` kwarg, not for any precision requirement
- Con: scientific notation (`%1.9e`) is avoided by IQ-TREE explicitly because "some software does not handle number format like '1.234e-6'" (see `mtree.cpp:462`)

#### C. Fixed high precision (e.g. 8 decimal, matching Augur)

Hardcode a generous fixed-decimal default like `%1.8f`.

- Pro: compact, predictable column widths in output files
- Pro: matches the Nextstrain ecosystem (Augur uses `%1.8f`)
- Con: still lossy -- small branch lengths (1e-9, common in large SARS-CoV-2 trees with ~30k sites) truncate to zero with fixed decimal
- Con: still an arbitrary choice requiring justification
- Con: fixed-decimal is the wrong unit -- significant digits preserve relative precision across magnitudes, decimal digits do not

#### D. Keep 3 significant digits (status quo)

- Pro: none identified
- Con: worst precision of any tool surveyed -- the next-lowest is Biopython at 5 decimal places
- Con: 0.4% relative error on typical branch lengths, destroys small branches entirely
- Con: contradicts `util-newick` behavior in the same codebase

Recommendation: A (full precision). Serialization should be lossless by default. Truncation is a consumer/display concern, not a writer concern. The `NwkWriteOptions` fields already exist for callers who want explicit truncation.

### Axis 2: scope of the `float_to_digits` fallback fix

The 3-sigfig default exists in two places: `fn format_weight()` in `treetime-io` (`.or(Some(3))`) and `fn float_to_digits()` in `treetime-utils` (line 82 fallback). Fix one or both?

#### A. Fix both

Remove the fallback in both locations: the `.or(Some(3))` in `format_weight` and the `max_significant_digits(3)` default in `float_to_digits`.

- Pro: consistent behavior -- `None` means "no limit" everywhere, matching the API contract implied by `Option<u8>`
- Pro: the `float_to_digits` fallback is a contract violation: the caller passes `None` ("no constraint") but the function secretly applies a constraint. This is a bug regardless of the precision debate
- Pro: zero blast radius -- every non-test caller of `float_to_digits` / `float_to_significant_digits` already passes an explicit digit count (verified: `run_loop.rs`, `clock_model.rs`, `rtt_chart_render.rs`, `console.rs`, `console_metrics.rs`, `console_tables.rs`, `domain_agreement.rs`)
- Con: changes `float_to_digits(x, None, None)` semantics for any future caller that relies on the implicit 3

#### B. Fix only `format_weight`, leave `float_to_digits` fallback

- Pro: narrower change surface
- Con: leaves a broken API contract in the general-purpose formatter -- `None` secretly means `Some(3)`
- Con: any future caller of `float_to_digits(x, None, None)` silently gets truncated to 3 significant digits
- Con: inconsistency between the two layers -- `format_weight` says "no limit" but `float_to_digits` disagrees

Recommendation: A (fix both). The `float_to_digits` fallback is independently wrong as an API contract. `Option::None` should mean "no constraint". No existing caller relies on it.

### Axis 3: CLI configurability

Should precision be exposed as a CLI flag?

#### A. No CLI flag (default full precision, API-only override via `NwkWriteOptions`)

- Pro: simpler CLI surface
- Pro: full precision is correct for all scientific use cases; truncation is a display/post-processing concern
- Con: users wanting compact output for visual inspection have no recourse except post-processing

#### B. Add `--precision N` flag (like RAxML-NG `--precision`)

- Pro: 3 of 10 surveyed tools expose this (RAxML-NG, MrBayes, BEAST2)
- Pro: useful for interop with parsers that struggle with long decimals
- Con: adds a flag that rarely matters; the correct default should not need overriding
- Con: requires plumbing through `OutputArgs` -> `NwkWriteOptions`

#### C. Defer CLI flag to a separate ticket

- Pro: keeps this change focused on fixing the bad default
- Pro: `NwkWriteOptions` fields already exist for programmatic use; CLI wiring is independent work
- Pro: no user has requested configurable precision -- solve the concrete problem first
- Con: if a user does need truncated output, they must wait for a follow-up

Recommendation: C (defer). The fix is removing a bad default. CLI exposure is orthogonal and should be tracked separately if a use case arises.

### Combined recommendation

A1 + A2 + C3: full precision default, fix both fallback sites, defer CLI flag.

Two-line code change: remove `.or(Some(3))` in `nwk.rs:301`, remove the `max_significant_digits(3)` fallback in `float.rs:82`. Zero blast radius on existing callers. Aligns `treetime-io` with `util-newick`'s existing full-precision behavior.

## Follow-up: Graphviz DOT writer

The Graphviz DOT writer [[src](https://github.com/neherlab/treetime/blob/b6c643bd/packages/treetime-io/src/graphviz.rs#L189)] also calls `format_weight` with `NwkWriteOptions::default()`. DOT output is for visualization, not data exchange -- full-precision labels like `0.00123456789012345` on every edge make graphs unreadable. After the default changes to full precision, the DOT writer should use its own truncated precision (e.g. 6 significant digits):

```rust
format_weight(weight, &NwkWriteOptions {
  weight_significant_digits: Some(6),
  ..NwkWriteOptions::default()
})
```

## Impact

- All analysis commands (ancestral, timetree, clock, optimize, prune, mugration) gain full-precision Newick output
- Existing `NwkWriteOptions` API is unchanged -- callers setting explicit digit limits see no behavior change
- Output files grow slightly (a few extra characters per branch length)
- Downstream tools consuming v1 output receive accurate branch lengths without silent truncation

## Validation

1. Round-trip property test: parse a tree with known branch lengths, write it, re-parse, assert values match within 1 ulp
2. Golden master comparison: v1 output branch lengths match v0 output within v0's formatting precision (7 decimal places)
3. Verify `util-newick` standalone writer behavior is unchanged (already full-precision)
4. Update existing `float_to_digits` test: `float.rs:126` test case `default_behavior` asserts `float_to_digits(1.23456, None, None)` produces `"1.23"` -- must change to expect `"1.23456"` after removing the 3-sigfig fallback

## Related

- Source issue: [kb/issues/N-io-nwk-writer-3-sigfig-default-truncates-precision.md](../issues/N-io-nwk-writer-3-sigfig-default-truncates-precision.md)
- Newick annotation dialects report: [kb/reports/newick-annotation-dialects.md](../reports/newick-annotation-dialects.md)
- Output format selection (implemented): [kb/proposals/output-format-selection.md](output-format-selection.md)
