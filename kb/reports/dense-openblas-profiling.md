# Dense mode performance investigation

Date: 2026-04-30
Commit: `554f3ea8` (`docs(known-issues): update dense test regression to reflect fixed tests`)

## 1. Problem statement

Dense marginal ancestral reconstruction is known to be slow. Three open known-issues document aspects of this:

- [M-timetree-marginal-dense-mpox-slow](../issues/M-timetree-marginal-dense-mpox-slow.md) -- 250x wall-time slowdown on mpox (20 tips, 44 min vs 7 sec for flu). Suspected cause: convolution grid spacing blow-up from low clock rate.
- [N-optimize-dense-iteration-slow](../issues/N-optimize-dense-iteration-slow.md) -- each marginal pass costs O(N _ s^2 _ L) over the full alignment length L. v0 mitigates with column compression; v1 dense does not compress.
- Dense Fitch compression gap (RESOLVED) -- dense now uses `PartitionFitch::compress` -> `infer_gtr` -> `into_marginal_dense`, and resolves indels via shared `fitch_indel` logic.

The goal of this investigation is to profile dense ancestral reconstruction, identify where time is spent, and find actionable improvements.

## 2. Setup

### 2.1 Dataset

All experiments use the same dataset and command:

- Dataset: `flu/h3n2/500` (500 tips, ~1800 nt alignment)
- Command: `treetime ancestral --method-anc=marginal --dense=true --tree=data/flu/h3n2/500/tree.nwk --outdir=tmp/profile/ancestral data/flu/h3n2/500/aln.fasta.xz`

### 2.2 Build configuration

Two Cargo profiles used:

- `profiling` -- inherits from `release` (LTO, opt-level 3) but keeps `debug = 1` and `strip = false`, so `perf` can resolve symbols. Used for profiling runs.
- `release` -- production build. Used for wall-time benchmarks.

The `bp` / `build-profiling` command in [`dev/dev`](../../dev/dev) builds the profiling profile via Docker (where OpenBLAS is available as a prebuilt static library):

```bash
./dev/docker/run ./dev/dev bp treetime   # profiling binary
./dev/docker/run ./dev/dev br treetime   # release binary
```

The profiling binary is at [`.build/docker/profiling/treetime`](../../.build/docker/profiling/treetime). The release binary is copied to [`.out/treetime`](../../.out/treetime).

### 2.3 OpenBLAS source

The project uses prebuilt static OpenBLAS from [binarylandia/build_openblas](https://github.com/binarylandia/build_openblas), installed via [`dev/docker/files/install-openblas`](../../dev/docker/files/install-openblas). The [2024-11-07 release](https://github.com/binarylandia/build_openblas/releases/tag/2024-11-07_22-51-18) provides both variants for each target triplet:

- `openblas-static-threads-0.3.28` -- threaded (includes internal thread pool)
- `openblas-static-0.3.28` -- sequential (no thread pool)

### 2.4 Host kernel settings

`perf_event_paranoid = -1`, `kptr_restrict = 0` (required for full symbol resolution and kernel frame capture).

## 3. Experiments

### 3.1 Baseline profile (threaded OpenBLAS)

**Method.** Built the profiling binary linked against `openblas-static-threads-0.3.28` (the variant originally configured in the project). Collected a CPU profile with `perf record`, then extracted flat self-time with `perf report`. Separately, built a release binary and measured wall time with `hyperfine`.

```bash
# Profile (profiling build, perf)
./dev/docker/run ./dev/dev bp treetime
perf record -g --call-graph dwarf -F 999 \
  -o tmp/perf-dense-ancestral.data \
  -- .build/docker/profiling/treetime \
  ancestral --method-anc=marginal --dense=true \
  --tree=data/flu/h3n2/500/tree.nwk \
  --outdir=tmp/profile/ancestral \
  data/flu/h3n2/500/aln.fasta.xz

perf report -i tmp/perf-dense-ancestral.data \
  --stdio --no-children --sort=sym --no-call-graph

# Wall time (release build, hyperfine)
./dev/docker/run ./dev/dev br treetime
cp .out/treetime tmp/treetime-threads-release
hyperfine --warmup 1 --runs 5 \
  "tmp/treetime-threads-release ancestral --method-anc=marginal --dense=true \
  --tree=data/flu/h3n2/500/tree.nwk --outdir=tmp/profile/ancestral \
  data/flu/h3n2/500/aln.fasta.xz"
```

**Profile results.** ~5000 samples across all threads. Self-time, flat:

| Self % | Symbol                            | Category            |
| -----: | --------------------------------- | ------------------- |
|  49.5% | `blas_thread_server`              | BLAS idle spin      |
|   5.0% | `crossbeam_deque::Stealer::steal` | rayon work-stealing |
|   4.2% | `__ieee754_log_fma`               | libm log            |
|   3.0% | `argmax_first`                    | treetime-utils      |
|   2.0% | `softmax_with_log_norm`           | treetime-utils      |
|   1.7% | `get_branch_mutation_matrix`      | GTR inference       |
|   1.4% | `accumulate_mutation_counts`      | GTR inference       |
|   1.2% | `normalize_inplace`               | marginal dense      |
|   1.1% | `Alphabet::seq2prof`              | sequence init       |
|   1.1% | `rayon wait_until_cold`           | rayon idle          |
|   0.9% | `__ieee754_exp_fma`               | libm exp            |
|   0.8% | `dgemm_kernel_HASWELL`            | BLAS matmul         |

**Wall time.**

| Metric   | Value         |
| -------- | ------------- |
| Wall     | 811 ms +/- 23 |
| User CPU | 2075 ms       |

### 3.2 Sequential OpenBLAS

**Method.** Changed [`dev/docker/files/install-openblas`](../../dev/docker/files/install-openblas) to download `openblas-static-0.3.28` instead of `openblas-static-threads-0.3.28` (both from the same [binarylandia release](https://github.com/binarylandia/build_openblas/releases/tag/2024-11-07_22-51-18)). Rebuilt Docker image. Clean-rebuilt both profiling and release binaries. Verified `blas_thread_server` symbol is absent from the binary. Ran the profile and benchmark procedure from 3.1.

```bash
# Verify the symbol is gone
nm .build/docker/profiling/treetime | grep blas_thread_server
# (no output)

# Profile (same perf command as 3.1)
perf record -g --call-graph dwarf -F 999 \
  -o tmp/perf-dense-ancestral-seq.data \
  -- .build/docker/profiling/treetime \
  ancestral --method-anc=marginal --dense=true \
  --tree=data/flu/h3n2/500/tree.nwk \
  --outdir=tmp/profile/ancestral \
  data/flu/h3n2/500/aln.fasta.xz

perf report -i tmp/perf-dense-ancestral-seq.data \
  --stdio --no-children --sort=sym --no-call-graph

# Wall time (A/B comparison of both release binaries)
cp .out/treetime tmp/treetime-seq-release
hyperfine --warmup 1 --runs 5 \
  -n "sequential" "tmp/treetime-seq-release ancestral --method-anc=marginal --dense=true \
  --tree=data/flu/h3n2/500/tree.nwk --outdir=tmp/profile/ancestral \
  data/flu/h3n2/500/aln.fasta.xz" \
  -n "threaded" "tmp/treetime-threads-release ancestral --method-anc=marginal --dense=true \
  --tree=data/flu/h3n2/500/tree.nwk --outdir=tmp/profile/ancestral \
  data/flu/h3n2/500/aln.fasta.xz"

# Run test suite
./dev/docker/run ./dev/dev t
```

**Profile results.** ~4100 samples across all threads. Self-time, flat:

| Self % | Symbol                            | Category            |
| -----: | --------------------------------- | ------------------- |
|   9.8% | `crossbeam_deque::Stealer::steal` | rayon work-stealing |
|   8.7% | `__ieee754_log_fma`               | libm log            |
|   6.7% | `argmax_first`                    | treetime-utils      |
|   3.4% | `get_branch_mutation_matrix`      | GTR inference       |
|   2.9% | `softmax_with_log_norm`           | treetime-utils      |
|   2.7% | `normalize_inplace`               | marginal dense      |
|   2.6% | `accumulate_mutation_counts`      | GTR inference       |
|   2.3% | `rayon wait_until_cold`           | rayon idle          |
|   2.1% | `Alphabet::seq2prof`              | sequence init       |
|   1.8% | `ndarray Zip::map_collect_owned`  | ndarray             |
|   1.8% | `__ieee754_exp_fma`               | libm exp            |
|   1.7% | `dgemm_oncopy_HASWELL`            | BLAS pack           |
|   1.3% | `dgemm_kernel_HASWELL`            | BLAS matmul         |

**Wall time.**

| Metric   | Value         |
| -------- | ------------- |
| Wall     | 781 ms +/- 22 |
| User CPU | 868 ms        |

**Test verification.** 2779 tests pass. 1 pre-existing known failure (`test_optimize_contribution_dense_sparse_ambiguous_r_value_and_gradient_consistency`) is unrelated to this change.

## 4. Conclusions

### 4.1 Threaded OpenBLAS was a configuration defect

The baseline profile (3.1) shows 49.5% of all CPU cycles spent in `blas_thread_server` -- the OpenBLAS internal thread pool busy-waiting for work. Actual BLAS compute (`dgemm_kernel_HASWELL`) accounts for 0.8%. Treetime's BLAS usage involves 5x5 matrices (nucleotide transition probability matrices). Threading overhead for matrices this small is pure waste.

Switching to sequential OpenBLAS (3.2) eliminates `blas_thread_server` entirely. Comparing the two runs:

| Metric    | Threaded (3.1) | Sequential (3.2) | Change            |
| --------- | -------------- | ---------------- | ----------------- |
| Wall time | 811 ms +/- 23  | 781 ms +/- 22    | 1.04x faster (4%) |
| User CPU  | 2075 ms        | 868 ms           | 2.4x fewer cycles |

Wall time improvement is modest because rayon already parallelizes real work across cores -- the BLAS threads were wasting cycles, not blocking the critical path. The CPU savings are the real win: 2.4x reduction in total user CPU, meaning less power, less thermal pressure, and less interference with other processes on the machine.

### 4.2 Where time goes after the fix

With the BLAS spin eliminated, the profile from 3.2 reveals where dense marginal reconstruction actually spends time. Filtering to application and library functions, normalized to 100% of useful work:

| Category                       | Share | Top functions                                                                            |
| ------------------------------ | ----: | ---------------------------------------------------------------------------------------- |
| Math primitives (log/exp/norm) | 36.2% | `log_fma`, `softmax_with_log_norm`, `normalize_inplace`, `exp_fma`, `normalize_from_log` |
| GTR inference                  | 13.0% | `get_branch_mutation_matrix`, `accumulate_mutation_counts`                               |
| Argmax / sequence extract      | 12.6% | `argmax_first`                                                                           |
| ndarray ops                    |  9.3% | `zip_mut_with`, `map_collect_owned`, `map`, `sum_axis`                                   |
| BLAS matmul                    |  6.5% | `dgemm_kernel`, `dgemm_oncopy`                                                           |
| Sequence init                  |  4.5% | `Alphabet::seq2prof`                                                                     |
| Map lookups / hashing          |  3.6% | `IndexMap::get`, `hash_one`                                                              |
| Other treetime                 |  3.3% | `assign_sequence`, `edge_subs`                                                           |

The dominant cost is log-space probability arithmetic in the marginal forward/backward passes: `log` and `exp` calls for normalization, softmax, and profile updates. This is intrinsic to the dense algorithm which processes every alignment column at every node.

### 4.3 Dense performance on this dataset is acceptable

781 ms wall time for 500 tips with full marginal reconstruction is not a bottleneck. The known-issue slowdowns (mpox 250x, optimize iteration cost) stem from other causes (convolution grid blow-up, lack of column compression) not triggered by this dataset.

### 4.4 Remaining overhead is algorithmic

Log/exp arithmetic (36%), per-column GTR inference (13%), and per-column argmax (13%) are all O(L) in alignment length. The primary lever for further improvement is reducing L through column compression. Rayon work-stealing at 9.8% suggests task granularity could also be coarsened.

## 5. Out-of-Docker builds

The Docker fix (switching to `openblas-static-0.3.28`) only affects Docker-based builds. Developers building on the host are affected by a separate configuration path.

### 5.1 How host builds link OpenBLAS

The workspace Cargo configuration uses `openblas-src` with the `system` feature:

```toml
# Cargo.toml (workspace)
openblas-src = { version = "=0.10.12", features = ["system"] }
blas-src = { version = "=0.12.1", features = ["openblas-src"] }
ndarray-linalg = { version = "=0.17.0", features = ["openblas-system"] }
```

The `system` feature means "link to whatever OpenBLAS is installed on the system via `pkg-config`." The build does not compile OpenBLAS from source. Which variant gets linked depends entirely on which system package the developer installed.

### 5.2 What the developer guide recommends

The [developer guide](../../docs/dev/developer_guide.md) instructs developers to install `libopenblas-dev`:

```bash
sudo apt-get install libopenblas-dev
```

On Ubuntu, `libopenblas-dev` is a metapackage with this dependency:

```
Depends: libopenblas-pthread-dev | libopenblas-openmp-dev | libopenblas-serial-dev
```

`apt` resolves alternatives left-to-right, so developers get `libopenblas-pthread-dev` -- the **threaded** variant. This means every developer following the guide gets the same threaded OpenBLAS that caused the 49.5% CPU waste observed in experiment 3.1.

### 5.3 Checking the current setup

**Check which system OpenBLAS package is installed:**

```bash
dpkg -l | grep openblas
```

- `libopenblas-pthread-dev` -- threaded (affected by this issue)
- `libopenblas-serial-dev` -- sequential (not affected)
- `libopenblas-openmp-dev` -- OpenMP threaded (also affected)

**Check if `OPENBLAS_LIB_DIR` overrides the system package:**

```bash
echo "${OPENBLAS_LIB_DIR:-not set}"
```

If set, the build uses the library at that path instead of the system package.

**Check a built binary for the threaded BLAS symbol:**

```bash
nm target/release/treetime | grep blas_thread_server
```

- Output with an address (e.g. `000... t blas_thread_server`) -- threaded OpenBLAS is linked. The binary is affected.
- No output -- sequential OpenBLAS is linked. The binary is not affected.

**Check the runtime environment variable:**

```bash
echo "${OPENBLAS_NUM_THREADS:-not set}"
```

If set to `1`, threaded OpenBLAS is capped to a single thread at runtime, mitigating the issue even if the threaded variant is linked.

### 5.4 Mitigation options

Three approaches, from least to most invasive:

1. **Runtime environment variable.** `OPENBLAS_NUM_THREADS=1` disables OpenBLAS internal threading regardless of which variant is linked. Developers can set it in their shell profile, or the treetime binary can set it at startup via `std::env::set_var` before any BLAS call. This is a safety net that works for all users on all platforms without any installation changes.

2. **Developer guide change.** Replace `libopenblas-dev` with `libopenblas-serial-dev` in the [developer guide](../../docs/dev/developer_guide.md). This installs the sequential variant system-wide. Developers who already have `libopenblas-dev` installed can switch with:

   ```bash
   sudo apt-get install libopenblas-serial-dev
   ```

   This only helps new setups and developers who re-read the guide.

3. **Use the binarylandia static library on host.** The same sequential `openblas-static-0.3.28` archive used by Docker can be downloaded and used for host builds. The `openblas-src` crate with `features = ["system"]` respects the `OPENBLAS_LIB_DIR` environment variable (set in [native.dockerfile:76](../../dev/docker/native.dockerfile#L76)). A developer can bypass the system package entirely:

   ```bash
   mkdir -p ~/.local
   curl -fsSL "https://github.com/binarylandia/build_openblas/releases/download/\
   2024-11-07_22-51-18/openblas-static-0.3.28-x86_64-unknown-linux-gnu-\
   2024-11-07_22-51-18.tar.xz" | tar -xJ -C ~/.local

   # Add to shell profile
   export OPENBLAS_LIB_DIR="$HOME/.local/lib"
   ```

   This statically links the sequential variant, matching the Docker build exactly. It also removes the system `libopenblas-dev` dependency.

### 5.5 Verifying a mitigation

After applying any mitigation from 5.4, rebuild from clean and re-check:

```bash
cargo clean
cargo build --release --bin=treetime
nm target/release/treetime | grep blas_thread_server
```

For mitigations 2 and 3 (system package change or binarylandia library), `blas_thread_server` should be absent from `nm` output. For mitigation 1 (runtime environment variable), the symbol may still be present, but the thread pool is disabled at runtime. To confirm, profile a short run:

```bash
perf record -F 99 -- target/release/treetime ancestral \
  --method-anc=marginal --dense=true \
  --tree=data/flu/h3n2/20/tree.nwk \
  --outdir=tmp/verify \
  data/flu/h3n2/20/aln.fasta.xz
perf report --stdio --no-children --sort=sym --no-call-graph | grep blas_thread_server
```

If `blas_thread_server` does not appear in the profile output, the mitigation is effective.

## 6. Proposals

### 6.1 Immediate (apply now)

- **Switch Docker builds to sequential OpenBLAS** -- change [`dev/docker/files/install-openblas`](../../dev/docker/files/install-openblas) to download `openblas-static-0.3.28` instead of `openblas-static-threads-0.3.28`. No code changes required. Eliminates 49.5% CPU waste for Docker-based builds.
- **Set `OPENBLAS_NUM_THREADS=1` at startup** -- add a one-line `std::env::set_var("OPENBLAS_NUM_THREADS", "1")` early in `main()`. Mitigates threaded OpenBLAS for all users on all platforms, including host builds and pre-built binaries distributed to end users.
- **Update developer guide** -- replace `libopenblas-dev` with `libopenblas-serial-dev` in the [developer guide](../../docs/dev/developer_guide.md) dependency list.

### 6.2 Near-term investigation

- **Column compression for dense mode** -- group identical alignment columns and weight by multiplicity. This is what v0 does. Would reduce L (and therefore wall time of marginal passes, GTR inference, and argmax) proportionally to alignment redundancy. For viral datasets, unique column count is typically 5-20% of L.
- **Profile the mpox dataset** -- the 250x slowdown documented in [M-timetree-marginal-dense-mpox-slow](../issues/M-timetree-marginal-dense-mpox-slow.md) is suspected to come from convolution grid blow-up in the timetree backward pass, not from ancestral reconstruction. A separate profiling session targeting `timetree --dense=true` on mpox would confirm or refute this.

### 6.3 Longer-term

- **Fast log/exp approximations** -- libraries like SLEEF or hand-tuned polynomial approximations could reduce the 36% math primitive cost, at some precision trade-off. Needs error analysis for the marginal inference use case.
- **Rayon task granularity** -- 9.8% in work-stealing suggests tasks are too fine-grained. Batching nodes or columns into larger work units could reduce synchronization overhead.
