# Performance: head vs base (before the refactor)

Date: 2026-10-07

The current code (`head`) is faster than the code before the refactor (`base`) on all 14 benchmark cases, and it uses less memory on 12 of them.

- **Time**: wall time, the time the user waits for a run, median of 5 runs, all threads
- **Memory**: peak RSS (maximum resident memory of the process), median of 5 runs
- **`base`**: commit `b39d2735` (2026-09-09), the last commit before the refactor
- **`head`**: `rust` tip `136647d1` (2026-10-07)

## Results

Negative change = `head` is better. Each table is sorted by the `base` value, largest first.

### Time

| Case | `base` | `head` | Change |
|---|--:|--:|--:|
| optimize dense `rsv/a/2000` | 40.5 s | 26.6 s | **-34%** |
| timetree `dengue/2000` | 35.3 s | 34.0 s | **-4%** |
| ancestral dense `rsv/a/2000` | 10.4 s | 3.1 s | **-71%** |
| timetree `rsv/a/2000` | 10.0 s | 9.3 s | **-7%** |
| optimize dense `dengue/500` | 5.9 s | 3.9 s | **-34%** |
| timetree `sc2/4500` | 4.8 s | 3.5 s | **-29%** |
| timetree `mpox/clade-ii/2000` | 3.7 s | 2.3 s | **-37%** |
| mugration `dengue/2000` | 2.9 s | 2.7 s | **-8%** |
| ancestral sparse `sc2/4500` | 1.54 s | 0.92 s | **-40%** |
| optimize sparse `dengue/2000` | 1.06 s | 0.76 s | **-28%** |
| timetree `mpox/clade-ii/500` | 0.83 s | 0.53 s | **-36%** |
| ancestral sparse `dengue/2000` | 0.35 s | 0.24 s | **-31%** |
| prune `dengue/2000` | 0.19 s | 0.11 s | **-42%** |
| clock `dengue/2000` | 0.06 s | 0.05 s | **-17%** |

### Memory

| Case | `base` | `head` | Change |
|---|--:|--:|--:|
| ancestral dense `rsv/a/2000` | 8002 MB | 7781 MB | **-3%** |
| optimize dense `rsv/a/2000` | 7985 MB | 7731 MB | **-3%** |
| timetree `mpox/clade-ii/2000` | 2308 MB | 1258 MB | **-46%** |
| timetree `dengue/2000` | 1784 MB | 1740 MB | **-2%** |
| timetree `sc2/4500` | 1658 MB | 1343 MB | **-19%** |
| optimize dense `dengue/500` | 1422 MB | 1293 MB | **-9%** |
| timetree `rsv/a/2000` | 1170 MB | 1060 MB | **-9%** |
| ancestral sparse `sc2/4500` | 1145 MB | 919 MB | **-20%** |
| timetree `mpox/clade-ii/500` | 568 MB | 324 MB | **-43%** |
| ancestral sparse `dengue/2000` | 295 MB | 245 MB | **-17%** |
| optimize sparse `dengue/2000` | 289 MB | 248 MB | **-14%** |
| prune `dengue/2000` | 118 MB | 79 MB | **-33%** |
| $\color{red}{\text{mugration dengue/2000}}$ | $\color{red}{\text{53 MB}}$ | $\color{red}{\text{96 MB}}$ | $\color{red}{\textbf{+81\%}}$ |
| $\color{red}{\text{clock dengue/2000}}$ | $\color{red}{\text{27 MB}}$ | $\color{red}{\text{28 MB}}$ | $\color{red}{\textbf{+4\%}}$ |

## Findings

- **Largest gains**: dense ancestral (-71% time), dense optimize (-34%), and timetree on `mpox` (-37% time, -46% memory)
- **Smallest gains**: the long timetree runs on `dengue/2000` (-4%) and `rsv/a/2000` (-7%). Most of their time goes to time inference, which changed little
- **Memory**: `head` is at or below `base` on every case except mugration and clock
- **Mugration memory grew** from 53 MB to 96 MB. The amount is small, but it is the only real memory regression. The cause is not investigated
- **Clock**: its changes are measurement noise (runs of 0.05 s and 27 MB)

## Output equivalence

Some results differ between `base` and `head`, so this comparison is close but not exact:

- Timetree values differ slightly (clock rate about 1e-10 relative on `mpox/500`, 1e-4 on `dengue/2000`)
- Dense optimize on `rsv/a/2000` writes different trees
- Ancestral now also reports mutations to ambiguous characters (e.g. `C7536M` on `sc2/4500`), so `head` writes more output there
- Sparse optimize, mugration, prune, and clock give the same results

## Details

### Builds

Each binary is the shipped Linux build of its commit: `CROSS_COMPILE=x86_64-unknown-linux-gnu ./dev/docker/run dev/cross/build treetime`, with `-C target-cpu=haswell`. That is the `release` profile for `base` and the `dist` profile for `head`.

### Method

- One idle x86-64 Linux machine, with no other load during the runs
- 5 rounds, binary order alternating every round. All runs exited with code 0
- `head` was measured in two separate series; they agree within 1% on time and 3% on memory, so the numbers are stable
- Each run was measured with `/usr/bin/time`, which also records CPU time (user + system); CPU time is not shown here
- Instruction counts were not collected

### Cases

Every case also has `--output-all=<dir>`. Timetree cases also have `--output-selection=nwk,nexus,auspice,augur-node-data,clock-model,gtr`.

| Case                           | Command                                                                                                                                                                     |
| ------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| optimize dense `rsv/a/2000`    | `optimize --tree=data/rsv/a/2000/tree.nwk --aln=data/rsv/a/2000/aln.fasta.xz --dense=true`                                                                                  |
| timetree `dengue/2000`         | `timetree --tree=data/dengue/2000/tree.nwk --dates=data/dengue/2000/metadata.tsv --aln=data/dengue/2000/aln.fasta.xz --name-column=genbank_accession --seed=1`              |
| ancestral dense `rsv/a/2000`   | `ancestral --tree=data/rsv/a/2000/tree.nwk --aln=data/rsv/a/2000/aln.fasta.xz --method-anc=marginal --dense=true --model=jc69`                                              |
| timetree `rsv/a/2000`          | `timetree --tree=data/rsv/a/2000/tree.nwk --dates=data/rsv/a/2000/metadata.tsv --aln=data/rsv/a/2000/aln.fasta.xz --name-column=accession --seed=1`                         |
| optimize dense `dengue/500`    | `optimize --tree=data/dengue/500/tree.nwk --aln=data/dengue/500/aln.fasta.xz --dense=true`                                                                                  |
| timetree `sc2/4500`            | `timetree --tree=data/sc2/4500/tree.nwk --dates=data/sc2/4500/metadata.tsv.xz --aln=data/sc2/4500/aln.fasta.xz --name-column=strain --seed=1`                               |
| timetree `mpox/clade-ii/2000`  | `timetree --tree=data/mpox/clade-ii/2000/tree.nwk --dates=data/mpox/clade-ii/2000/metadata.tsv --aln=data/mpox/clade-ii/2000/aln.fasta.xz --name-column=accession --seed=1` |
| mugration `dengue/2000`        | `mugration --tree=data/dengue/2000/tree.nwk --metadata=data/dengue/2000/metadata.tsv --attribute=country --metadata-id-columns=genbank_accession`                           |
| ancestral sparse `sc2/4500`    | `ancestral --tree=data/sc2/4500/tree.nwk --aln=data/sc2/4500/aln.fasta.xz --method-anc=marginal --dense=false`                                                              |
| optimize sparse `dengue/2000`  | `optimize --tree=data/dengue/2000/tree.nwk --aln=data/dengue/2000/aln.fasta.xz --dense=false`                                                                               |
| timetree `mpox/clade-ii/500`   | `timetree --tree=data/mpox/clade-ii/500/tree.nwk --dates=data/mpox/clade-ii/500/metadata.tsv --aln=data/mpox/clade-ii/500/aln.fasta.xz --name-column=accession --seed=1`    |
| ancestral sparse `dengue/2000` | `ancestral --tree=data/dengue/2000/tree.nwk --aln=data/dengue/2000/aln.fasta.xz --method-anc=marginal --dense=false`                                                        |
| prune `dengue/2000`            | `prune --tree=data/dengue/2000/tree.nwk --aln=data/dengue/2000/aln.fasta.xz --prune-short=1e-6 --prune-empty`                                                               |
| clock `dengue/2000`            | `clock --tree=data/dengue/2000/tree.nwk --dates=data/dengue/2000/metadata.tsv --name-column=genbank_accession`                                                              |

### Reproduce

Build each commit as above, then run each case with both binaries from the repository root and wrap every run in `/usr/bin/time`, for example:

```bash
/usr/bin/time -f '%e s %M KB' .out/treetime-x86_64-unknown-linux-gnu timetree --tree=data/sc2/4500/tree.nwk --dates=data/sc2/4500/metadata.tsv.xz --aln=data/sc2/4500/aln.fasta.xz --name-column=strain --seed=1 --output-selection=nwk,nexus,auspice,augur-node-data,clock-model,gtr --output-all=tmp/bench/timetree-sc2-4500
```

## Open items

- Profile the mugration memory growth
- Collect instruction counts with `perf stat`, which do not depend on machine load
