# Newick date annotation uses 2-decimal precision

v1 formats the `date` field in Newick/Nexus node comments with 2 decimal places (`{time:.2}`), matching v0's `%1.2f` format. Two decimal places gives ~3.65 day resolution (0.01 year), which is coarse for fast-evolving pathogens sampled days apart.

## v0 behavior

v0 uses different precision for the same `numdate` value depending on output format:

- Newick/Nexus comment: `%1.2f` (2 decimals) at `CLI_io.py`, line `n.comment += ... + 'date=%1.2f' % n.numdate`
- TSV dates file: `%f` (6 decimals, Python default) at `CLI_io.py`, line `fh_dates.write('%s\t%s\t%f\n' % (n.name, n.date, n.numdate))`
- Auspice JSON: full `float()` precision at `CLI_io.py`, line `j['node_attrs']['num_date'] = {'value': float(n.numdate)}`

The Newick format is the lowest-fidelity output.

## v1 behavior

v1 writes the annotations of a node in the same text as v0, for example `[&mutations="A55G,T93C",date=2003.84]`:

- **Key order**: `mutations`, then `date`, then the trait attribute of `mugration`, the order in which v0 builds the comment
- **Quoting**: every string value is written in double quotes in BEAST style, with an embedded `"` doubled, as v0 writes `&mutations="..."` and `&<attribute>="<value>"` ([packages/legacy/treetime/treetime/wrappers.py#L903](../../packages/legacy/treetime/treetime/wrappers.py#L903)). NHX values stay unquoted
- **Date**: two decimals, the text of v0's `%1.2f`, written unquoted

`fn nwk_node_comments()` in [packages/app-output/src/nwk_comments.rs](../../packages/app-output/src/nwk_comments.rs) builds this list with typed values, so the writer never guesses a value's type from its text: the date is `NewickValue::NumberText`, which keeps the trailing zero of `2020.50`, and trait values such as `01`, `Nan` or `true` stay strings.

## Rationale

The 2-decimal format is a deliberate v0 design choice for Newick readability. Newick strings are human-inspectable, and 2 decimals keeps annotations compact. Higher-fidelity date output is available in other formats (TSV, JSON). Matching v0 avoids downstream tool breakage for users parsing Newick date annotations with fixed-width expectations.

## Precision analysis

| Decimals | Resolution (days) | Resolution (hours) |
| -------- | ----------------- | ------------------ |
| 2        | 3.65              | 87.66              |
| 4        | 0.04              | 0.88               |
| 6        | 0.0004            | 0.01               |

## Impact

Newick/Nexus date annotations lose sub-week temporal resolution. Full-precision dates are available in the Auspice JSON output (`num_date` attribute). The TSV dates output is not yet implemented (see `kb/issues/N-timetree-node-dates-output-unimplemented.md`).
