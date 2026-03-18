# Timetree Test Matrix

Systematic testing of v1 `timetree` command across datasets and parameter combinations.

## Matrix

Rows are exact CLI args added to base command. Empty row = default (no extra args).

| Args                                          | flu/h3n2/20 | flu/h3n2/200 | flu/h3n2/500 | ebola/20 | zika/20 | rsv/a/20   | tb/20      | dengue/500 | lassa/L/20 | mpox/clade-ii/20 |
| --------------------------------------------- | ----------- | ------------ | ------------ | -------- | ------- | ---------- | ---------- | ---------- | ---------- | ---------------- |
| (default)                                     | OK          | OK           | OK           | OK       | OK      | CRASH-grid | CRASH-grid | CRASH-grid | CRASH-grid | CRASH-grid       |
| `--coalescent=10`                             | OK,coal     | OK,coal      | -            | OK,coal  | OK,coal | -          | -          | -          | -          | -                |
| `--coalescent=0.1`                            | OK,coal     | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--coalescent=1000`                           | OK,coal     | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--coalescent-skyline`                        | OK,coal     | OK,coal      | -            | OK,coal  | -       | -          | -          | -          | -          | -                |
| `--coalescent-opt`                            | OK          | OK           | -            | OK       | OK      | -          | -          | -          | -          | -                |
| `--time-marginal=always`                      | OK          | OK           | OK           | OK       | OK      | -          | -          | -          | -          | -                |
| `--coalescent-skyline --time-marginal=always` | OK,coal     | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--coalescent-opt --time-marginal=always`     | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--branch-length-mode=input`                  | OK          | -            | -            | OK       | -       | OK         | OK         | neg-err    | -          | -                |
| `--clock-rate=0.003`                          | OK          | -            | -            | OK       | -       | -          | -          | -          | -          | -                |
| `--max-iter=0`                                | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--max-iter=4`                                | OK          | OK           | -            | OK       | -       | -          | -          | -          | -          | -                |
| `--confidence`                                | OK,conf     | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--covariation`                               | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--reroot=min-dev`                            | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--reroot=oldest`                             | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--keep-root`                                 | CRASH-grid  | CRASH-grid   | -            | OK       | OK      | CRASH-grid | CRASH-grid | -          | -          | -                |
| `--gtr=jc69`                                  | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--keep-polytomies`                           | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--resolve-polytomies`                        | OK          | -            | -            | OK       | -       | -          | -          | -          | -          | -                |
| `--method-anc=parsimony`                      | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--confidence --clock-std-dev=0.0005`         | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--tip-slack=0.1`                             | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--relax=...`                                 | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |
| `--time-marginal=only-final`                  | OK          | -            | -            | -        | -       | -          | -          | -          | -          | -                |

## Legend

### Cell values

- `-` = not tested
- `OK` = passed
- `CRASH-<id>` = crash, see issue
- `HANG` = command timed out (deadlock, infinite loop)
- `<id>` or `<id>,<id>` = non-crash issues found

### Issue references

| ID      | Issue file                                                                                  |
| ------- | ------------------------------------------------------------------------------------------- |
| grid    | [H-timetree-crash-grid-zero-branch](H-timetree-crash-grid-zero-branch.md)                   |
| coal    | [M-timetree-coalescent-ci-excludes-internal](M-timetree-coalescent-ci-excludes-internal.md) |
| neg-err | Proper error on negative clock rate (data lacks temporal signal)                            |

## Base Command

```bash
timeout 120 ./dev/docker/run ./dev/dev r treetime -- timetree \
  --tree=data/$dataset/tree.nwk \
  --dates=data/$dataset/metadata.tsv \
  --outdir=tmp/timetree/$dataset \
  data/$dataset/aln.fasta.xz
```
