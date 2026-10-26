---
title: Testing against the reference
description: Compare the outputs of the CLI with a baseline snapshot, and run TreeTime v0 and the Python scripts.
---

`dev/smoke` (host, needs Docker and `zstd`) runs the CLI over a matrix of commands, datasets and flag variants, stores the outputs in a snapshot named after the source tree that built the binary, and compares them byte for byte with a baseline snapshot, by default the one of the tip of the base branch `rust`. It reports crashes, timeouts, missing outputs and changed outputs, with the command, stderr and the changed values of each case. It also lists the cases that pass in both snapshots but print a different stderr, such as a new or repeated warning; these differences are information and never fail a case. The cases are defined in `dev/smoke.toml`, one row per command variant with its flags, datasets and expected failures; `./dev/smoke --overview` prints the cases per command and variant. `./dev/smoke --help` describes every option, the snapshot layout and the exit codes.

- `just smoke` compares the quick tier with `rust`
- `just smoke-all` compares every case with `rust`
- `just smoke-run` runs the cases without a baseline: crash, timeout and output checks only
- `just smoke-failed` runs again the cases that did not pass in the last run
- `just smoke-prune` deletes the dirty snapshots of the checkout other than the current one, and old snapshots that no branch needs. Every run also deletes the earlier current-tree snapshots of its checkout (see [Snapshot store](#snapshot-store))

```bash
./dev/smoke                    # quick tier against rust
./dev/smoke --tier full        # every case
./dev/smoke --against <ref>    # another commit, branch or snapshot id as the baseline
./dev/smoke --no-compare       # statuses only
./dev/smoke --only <regex>     # cases whose id matches
./dev/smoke --rerun-failed     # cases that did not pass in the last run
./dev/smoke --list             # print the selected case ids
```

Other options:

- `--rerun`: run the selected cases again even when their stored results are reusable
- `--no-build`: use the existing `.out/treetime` for the current snapshot
- `--check-determinism`: run each case a second time with `--jobs=1` and mark the case nondeterministic when the outputs differ. Only the analysis commands and `pipeline` take `--jobs`; the cases of other commands run a second time unchanged
- `--jobs N`, `--timeout-scale X`: parallel cases, and a factor on every case timeout
- `--memory-budget GIB`: memory the running cases may use together, 75% of the available memory by default. Every treetime process is capped at the budget, so a case that needs more fails with an allocation error instead of exhausting the machine
- `--min-free GIB`, `--max-file-size GIB`: the disk limits (see [Disk limits](#disk-limits))

Tiers:

- `quick` (default): datasets of at most 100 sequences, the help texts and a few larger timetree cases, about 6 minutes
- `full`: every case, about 40 minutes

## Snapshot store

Without settings, each checkout keeps its own store: `snapshots/` at the checkout root (git-ignored), with the temporary baseline worktrees and the lock files under `tmp/smoke/`. With `SMOKE_STORE` set to an absolute directory in `.env` (see [Machine settings](/development/setup#machine-settings)), all checkouts on the machine share one store with `snapshots/`, `worktrees/` and `locks/` in that directory. A baseline is then built once per machine, and every worktree reuses its binary and case results.

- **Snapshot id**: the first 10 hex digits of the git tree id of the commit, or `<tree>+dirty-<hash>` for a working tree with uncommitted or untracked changes. The id depends on content only, so commits with the same tree, such as a branch and the target it was fast-forwarded into, share one snapshot. `--against` resolves a git ref to the tree of its commit, and also accepts an existing snapshot id
- **Binaries**: the current binary is built in the checkout unless its snapshot already has one. A git ref baseline is built in a temporary worktree of the store with a copy of the checkout's `.env`, so it uses the same compiler cache. `dev/smoke` also sets `TREETIME_CARGO_HOME` of `dev/docker/run` to the cargo home of the checkout, so the baseline build reuses the crates that the checkout downloaded; a baseline revision whose `dev/docker/run` predates the setting downloads them again
- **Reuse**: a stored case is reused when its command, expectation, declared outputs, working directory, input files, file size limit and output storage layout are unchanged. A snapshot whose binary changed loses its stored cases
- **Outputs**: an output file of at least 1 MiB is stored compressed as `<name>.zst`; read it with `zstd -dc`. `outputs.json` holds the SHA-256 of every output after removing the volatile fields, so a comparison decompresses a file only when its hash differs
- **Locks**: a run holds an exclusive lock on a snapshot while it builds and runs it, and a shared lock while it reads the results for its report. A second run that needs the same snapshot logs the pid and checkout of the holder, waits, and then reuses the stored results. Runs on different snapshots never wait for each other. The lock is released when its holder exits, also on a crash, and an interrupted build leaves no binary behind
- **Pruning after every run**: each edit or commit gives the current tree a new snapshot that one comparison uses. After every run, `dev/smoke` deletes the snapshots that the checkout ran before as its current tree, except the current one, the baseline of the run, the snapshots of local branch tips, and snapshots used as a baseline
- **Pruning on request**: `--prune` deletes the dirty snapshots that the current checkout wrote, other than its current one, and the snapshots of committed trees older than 7 days whose tree is not the tip of any local branch. Both kinds of pruning skip snapshots that another run holds, so one worktree never deletes the work of another

Layout of one snapshot:

```text
snapshots/<id>/
  snapshot.json             commit, tree, dirty hash, ref, checkout, binary SHA-256, creation time
  bin/treetime              the binary that produced the snapshot
  last-run.json             statuses of the last run and its checkout (for --rerun-failed)
  cases/<case-id>/
    out/                    outputs; a file of at least 1 MiB is stored compressed as <name>.zst
    outputs.json            SHA-256 of every output after removing the volatile fields, and whether it is compressed
    stdout.txt, stderr.txt
    cmd.sh                  the exact command, runnable from the checkout root in the container
    result.json             status, exit code, wall and CPU time, peak memory, case key
  compare/<baseline-id>/    report.md and report.tsv of a comparison
  report.md, report.tsv     report of a --no-compare run
```

`report.md` groups failing cases by likely cause and lists value-level differences. Before the comparison, only the volatile fields `meta.updated` of `*.auspice.json` and `generated_by.version` of `*.node-data.json` are removed; numeric differences are reported, never tolerated. A `changed` case, a `regressed` case, a `still-failing` case without a declared expectation, and an `unexpected-pass` case fail the run.

Exit codes: 0 when every selected case passed, 1 when a case failed or changed, 2 when the script could not complete (build failure, git error, invalid arguments, container failure, free disk space below the limit).

## Disk limits

A full run writes large outputs; the `.dot` files of the larger datasets reach about 1 GB each. Two limits keep a run from filling the disk:

- **Free disk space** (`--min-free GIB`, default `SMOKE_MIN_FREE_GIB` in `.env`, else 20): no build starts while the file system of the checkout or the store has less free space. While cases run, the free space of the store is checked every second; below the limit, the running cases are stopped, their partial outputs deleted, no further case starts, and the run exits with code 2. The cases that completed keep their results, so the next run continues from them
- **File size** (`--max-file-size GIB`, default 8): every treetime process may write files up to this size. A larger file stops the process with `SIGXFSZ`, and the case fails with that reason

## TreeTime v0 and Python scripts

TreeTime v0 (`packages/legacy/treetime`), the golden-master capture scripts, and the research scripts in `test_scripts/` share one Python environment, defined in `pixi.toml` and locked in `pixi.lock`. It pins every package and uses the OpenBLAS 0.3.34 release of conda-forge on a single thread. The Rust build links OpenBLAS develop `de5d4914`, 432 commits after 0.3.34, with LAPACK compiled by flang, so v0 and v1 run close but not identical linear algebra kernels. The environment puts v0 of the checkout on `PYTHONPATH`, so both the v0 `treetime` command and `import treetime` in scripts work.

In the Python image:

```bash
./dev/docker/python treetime ancestral --help
./dev/docker/python python3 test_scripts/fitch.py
```

On the host, [pixi](https://pixi.sh) comes with `mise install` and creates the environment in `.pixi/` on first use:

```bash
pixi run treetime ancestral --help
pixi run python3 test_scripts/fitch.py
```

After changing `pixi.toml`, run `pixi lock` and commit both files. Package releases must be at least seven days old; `exclude-newer` in `pixi.toml` enforces this.
