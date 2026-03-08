# Internal node dates missing with bad fixed clock rate

When `--clock-rate` is set to a value far from the estimated rate, internal nodes
lose date annotations. The solver produces degenerate time distributions when
the rate is wrong.

## Example

ebola/20 estimated rate: ~0.000626

- `--clock-rate=0.0006` (close): 39/39 annotations (all present)
- `--clock-rate=0.003` (5x too high): 20/39 annotations (tips only)

## Root cause

Same mechanism as
[Internal node dates missing at scale](M-timetree-internal-dates-missing-scale.md):
degenerate time distributions cause `likely_time()` to return `None`. A wrong
clock rate shifts the expected branch lengths, producing distributions that are
too narrow or too wide for the solver to find a valid mode.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 --clock-rate=0.003 \
  --tree=data/ebola/20/tree.nwk --dates=data/ebola/20/metadata.tsv \
  --outdir=tmp/repro-bad-rate data/ebola/20/aln.fasta.xz
grep -oP '\[&[^\]]*\]' tmp/repro-bad-rate/timetree.nexus | wc -l
# 20 (expected: 39)
```

## Related issues

- [Negative clock rate causes assertion failure](H-timetree-negative-clock-rate-assertion.md)
  extreme case where the rate is not just wrong but negative
