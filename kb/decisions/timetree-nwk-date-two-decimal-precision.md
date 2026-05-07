# Newick date annotation uses 2-decimal precision

v1 formats the `date` field in Newick/Nexus node comments with 2 decimal places (`{time:.2}`), matching v0's `%1.2f` format. Two decimal places gives ~3.65 day resolution (0.01 year), which is coarse for fast-evolving pathogens sampled days apart.

## v0 behavior

v0 uses different precision for the same `numdate` value depending on output format:

- Newick/Nexus comment: `%1.2f` (2 decimals) at `CLI_io.py`, line `n.comment += ... + 'date=%1.2f' % n.numdate`
- TSV dates file: `%f` (6 decimals, Python default) at `CLI_io.py`, line `fh_dates.write('%s\t%s\t%f\n' % (n.name, n.date, n.numdate))`
- Auspice JSON: full `float()` precision at `CLI_io.py`, line `j['node_attrs']['num_date'] = {'value': float(n.numdate)}`

The Newick format is the lowest-fidelity output.

## v1 behavior

v1 matches v0's Newick format exactly:

```rust
comments.insert("date".to_owned(), format!("{time:.2}"));
```

At `packages/treetime/src/representation/payload/timetree.rs#L135`.

## Rationale

The 2-decimal format is a deliberate v0 design choice for Newick readability. Newick strings are human-inspectable, and 2 decimals keeps annotations compact. Higher-fidelity date output is available in other formats (TSV, JSON). Matching v0 avoids downstream tool breakage for users parsing Newick date annotations with fixed-width expectations.

## Precision analysis

| Decimals | Resolution (days) | Resolution (hours) |
| -------- | ----------------- | ------------------ |
| 2        | 3.65              | 87.66              |
| 4        | 0.04              | 0.88               |
| 6        | 0.0004            | 0.01               |

## Impact

Newick/Nexus date annotations lose sub-week temporal resolution. Full-precision dates are available in the Auspice JSON output (`num_date` attribute). The TSV dates file (`write_node_dates`) is not yet implemented (see `kb/tickets/timetree-output-implement-node-dates.md`).
