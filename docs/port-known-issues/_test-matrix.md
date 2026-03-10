# Timetree Test Matrix

Systematic testing of v1 `timetree` command across datasets and parameter combinations.

## Matrix

Rows are exact CLI args added to base command. Empty row = default (no extra args).

| Args                                   | flu/h3n2/20    | flu/h3n2/200   | flu/h3n2/500     | ebola/20       | zika/20        | rsv/a/20   | tb/20      | dengue/500 | lassa/L/20 | mpox/clade-ii/20 |
| -------------------------------------- | -------------- | -------------- | ---------------- | -------------- | -------------- | ---------- | ---------- | ---------- | ---------- | ---------------- |
| (default)                              | OK\*           | OK\*           | OK\*,dates-scale | OK\*           | OK\*           | CRASH-grid | CRASH-grid | CRASH-grid | CRASH-grid | CRASH-grid       |
| `--coalescent=10`                      | OK\*,coal,rate | OK\*,coal,rate | coal,dates-scale | OK\*,coal,rate | OK\*,coal,rate | -          | -          | -          | -          | -                |
| `--coalescent=0.1`                     | OK\*,coal      | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--coalescent=1000`                    | OK\*,coal      | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--coalescent=skyline`                 | OK\*,coal      | OK\*,coal      | -                | OK\*,coal      | -              | -          | -          | -          | -          | -                |
| `--coalescent-opt`                     | OK\*           | OK\*,rate      | -                | OK\*,rate      | OK\*,rate      | -          | -          | -          | -          | -                |
| `--time-marginal`                      | OK\*           | OK\*           | OK\*,dates-scale | OK\*,rate      | OK\*,rate      | -          | -          | -          | -          | -                |
| `--time-marginal=always`               | OK\*,always    | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--coalescent=skyline --time-marginal` | OK\*,coal,ci   | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--coalescent-opt --time-marginal`     | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--branch-length-mode=input`           | OK\*,coal      | -              | -                | OK\*,coal      | -              | OK\*,neg   | OK\*,coal  | OK\*,neg   | -          | -                |
| `--clock-rate=0.003`                   | OK\*           | -              | -                | OK\*,dates-bad | -              | -          | -          | -          | -          | -                |
| `--max-iter=0`                         | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--max-iter=4`                         | OK\*           | OK\*,rate      | -                | OK\*,rate      | -              | -          | -          | -          | -          | -                |
| `--confidence`                         | OK\*,conf      | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--covariation`                        | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--reroot=min-dev`                     | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--reroot=oldest`                      | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--keep-root`                          | CRASH-grid     | CRASH-grid     | -                | OK\*           | OK\*           | CRASH-grid | CRASH-grid | -          | -          | -                |
| `--gtr=jc69`                           | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--keep-polytomies`                    | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--resolve-polytomies`                 | OK\*           | -              | -                | OK\*           | -              | -          | -          | -          | -          | -                |
| `--method-anc=parsimony`               | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--vary-rate`                          | CRASH-vary     | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--tip-slack=0.1`                      | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--relax=...`                          | CRASH-relax    | -              | -                | -              | -              | -          | -          | -          | -          | -                |
| `--only-final`                         | OK\*           | -              | -                | -              | -              | -          | -          | -          | -          | -                |

## Legend

### Cell values

- `-` = not tested
- `OK` = passed
- `OK*` = passed with `--clock-filter=0` workaround (bypasses deadlock)
- `CRASH-<id>` = crash, see issue
- `HANG` = command timed out (deadlock, infinite loop)
- `<id>` or `<id>,<id>` = non-crash issues found

### Issue references

| ID          | Issue file                                                                                  |
| ----------- | ------------------------------------------------------------------------------------------- |
| grid        | [H-timetree-marginal-dense-backward-crash](H-timetree-marginal-dense-backward-crash.md)     |
| vary        | [H-timetree-vary-rate-unimplemented](H-timetree-vary-rate-unimplemented.md)                 |
| relax       | [M-timetree-relax-arg-parsing](M-timetree-relax-arg-parsing.md)                             |
| coal        | [M-timetree-coalescent-ci-excludes-internal](M-timetree-coalescent-ci-excludes-internal.md) |
| rate        | (clock rate frozen across iterations - not yet documented)                                  |
| conf        | [M-timetree-confidence-flag-ignored](M-timetree-confidence-flag-ignored.md)                 |
| always      | [M-timetree-time-marginal-always-ignored](M-timetree-time-marginal-always-ignored.md)       |
| dates-scale | [M-timetree-internal-dates-missing-scale](M-timetree-internal-dates-missing-scale.md)       |
| dates-bad   | [M-timetree-internal-dates-bad-fixed-rate](M-timetree-internal-dates-bad-fixed-rate.md)     |
| ci          | [M-timetree-coalescent-ci-excludes-internal](M-timetree-coalescent-ci-excludes-internal.md) |
| neg         | (negative clock rate silently stored in input-BL mode - not yet documented)                 |

## Base Command

```bash
timeout 120 ./dev/docker/run ./dev/dev r treetime -- timetree \
  --tree=data/$dataset/tree.nwk \
  --dates=data/$dataset/metadata.tsv \
  --outdir=tmp/timetree/$dataset \
  data/$dataset/aln.fasta.xz
```

Add `--clock-filter=0` for `OK*` runs.

## Notes

- Rows with multiple args test specific interactions
- Add new rows as needed for additional parameter combinations
- Prefer testing on small datasets first, then expand to larger ones
