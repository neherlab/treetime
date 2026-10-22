# The outputs of timetree depend on the rounding of the platform

Differences in the last bits of linear algebra and math functions can grow through the time inference of `timetree` to relative differences in the clock rate above the project threshold of 1e-6. Changes of the OpenBLAS build have changed the clock rate by up to 1.7e-3 relative. Between the two OpenBLAS kernel sets of the current native build, the differences stay far below the threshold. `ancestral` gives byte-identical outputs on the CLI and in the WebAssembly back end of web2.

## Evidence

The current builds use OpenBLAS develop `de5d4914` on every target: natively with LAPACK compiled by flang 22, in web2 with the LAPACK C translation (the `openblas` entries of `.config/mise.toml` and `.config/mise.image.toml`).

Measured on the current code with this build, single-threaded, `timetree` with default settings and, for `data/zika/20`, the configuration that web2 resolves for a run with tree, alignment, and metadata (`--max-iter=2`, marginal branch lengths, `--clock-filter=3`, the default at the time of the measurement):

| Dataset | Output | CLI, OpenBLAS `Haswell` kernels | CLI, `OPENBLAS_CORETYPE=Prescott` (`Katmai` kernels) | Relative difference |
| --- | --- | --- | --- | --- |
| `zika/20` | `clock_rate` of `timetree.clock-model.json` | 0.0010255720768280348 | 0.001025572076832962 | 4.8e-12 |
| `zika/20` | `log_lh_time` of the first iteration (`tracelog.csv`) | -40.37081080491443 | -40.3708108065102 | 4.0e-11 |
| `zika/20` | `log_lh_seq` of the first iteration | -17281.96059624937 | -17281.96059625 | 3.6e-14 |
| `ebola/100` | `clock_rate` | 0.000728844005662043 | 0.0007288440068025652 | 1.6e-9 |
| `tb/20` | `clock_rate` | 0.0005375728626081336 | 0.0005375728625939621 | 2.6e-11 |

On `zika/20`, `timetree.nwk`, `timetree.nexus`, `timetree.auspice.json`, and `timetree.sequences.fasta` are byte-identical between the two kernel sets; `timetree.model.json`, `timetree.clock-model.json`, and `timetree.node-data.json` differ in the last digits. web2 is not measured with this build.

Earlier measurements show the amplification:

- With OpenBLAS 0.3.28 natively and OpenBLAS develop `eb00d597` in web2, `data/zika/20` with the web2 configuration gave clock rates 0.0010266268996163546 (`Haswell`), 0.0010266466684456434 (`Katmai`), and 0.001026646778721711 (web2, `WASM128_GENERIC`): 1.9e-5 relative between web2 and `Haswell`, and 1.1e-7 between web2 and `Katmai`. The log-likelihood of the sequences in the first iteration differed by 2.0e-7 relative between web2 and the CLI, and the root-to-tip divergence of the first tip by 1.2e-9 relative
- The OpenBLAS version alone changes `timetree` beyond the threshold. `dev/smoke` with the same code, once with OpenBLAS 0.3.28 (commit `2dabfc59`) and once with 0.3.34, both with `Haswell` kernels on an Intel i9-10900K: the clock rate changes by 9.2e-4 relative on `lassa/L/50`, 3.9e-4 on `ebola/362` with skyline, 3.2e-4 on `ebola/100`, 2.6e-4 on `tb/20`, and 7.9e-6 on `rsv/a/100`; the branch lengths of the time tree change by up to 0.4 on `lassa/L/50`, and the reconstructed sequences of `lassa/L/50` and of `ebola/362` with `--clock-filter=10` change. `timetree` on `zika/20` stays byte-identical. The outputs of `clock`, `homoplasy`, and `prune` stay byte-identical, and those of `ancestral`, `mugration`, and `optimize` change by at most 1.8e-7 relative, apart from the last printed digit of one value in `mugration.trait-probabilities.tsv`
- The OpenBLAS build changes `timetree` beyond the threshold again. `dev/smoke` (quick tier) of commit `2e275920` against `4f116c6d`, the same code once with the OpenBLAS 0.3.34 release of the earlier images, once with develop `de5d4914` and LAPACK compiled by flang 22, both with `Haswell` kernels on an Intel i9-10900K: `timetree` changes in 58 cases; the clock rate changes by 1.7e-3 relative on `ebola/20` with coalescent optimization, 1.5e-3 on `flu/h3n2/20` dense and with coalescent, 3.7e-4 on `dengue/20`, 3.4e-4 on `ebola/362` with skyline, 3.3e-4 on `ebola/100`, 3.2e-4 on `tb/20`, below 3e-4 on the other datasets, and by at most 2.4e-11 on most `flu/h3n2/20` variants. The other commands change by at most: `ancestral` 1.9e-12 relative (GTR of `--gtr-iterations=3`), `optimize` 4.3e-11 apart from one branch of 3.7e-8 on `ebola/100` (1.1e-5 relative, 4.2e-13 absolute), `mugration` 1.4e-7 apart from the last printed digit of one value in `mugration.trait-probabilities.tsv` and one tied state on `dengue/100`, which the tie rule of [kb/decisions/posterior-most-likely-state-tie-rule.md](../decisions/posterior-most-likely-state-tie-rule.md) now fixes. `clock`, `homoplasy`, and `prune` stay byte-identical. `cli/pipeline-steps-timetree` and `config/zika/20/pipeline` change through their `timetree` step
- `ancestral` on `data/flu/h3n2/20`: `node-data.json`, `model.json`, `nexus`, `nwk`, and `sequences.fasta` of web2 and of the CLI are byte-identical
- The compile-time CPU level of TreeTime's own code is not a rounding source: x86_64 builds with `-C target-cpu=haswell`, `x86-64-v3`, and `x86-64-v2` (liblzma and libzstd unchanged) wrote byte-identical outputs for the `flu/h3n2/500` pipeline, dense and sparse `ancestral`, `optimize`, `mugration`, and `timetree` on `lassa/L/500`, `rsv/a/2000`, and `flu/h3n2/500` with the skyline and confidence intervals, at 1 and at 20 threads. Rust does not fuse multiplications and additions into FMA without `mul_add`, which no TreeTime code calls
- Rounding sources that differ between the platforms: the OpenBLAS kernels (FMA in `Haswell`, none in `Katmai` and `WASM128_GENERIC`); LAPACK (OpenBLAS develop `de5d4914` compiled by flang 22 natively, its C translation with the LAPACK 3.9 `?lartg` and `?lassq` in web2); `exp`, `log`, and `pow` (glibc natively, with FMA variants on FMA CPUs; musl of Emscripten, from the same ARM optimized routines, in web2); the small products of `matrixmultiply` (FMA kernels natively, the fallback kernel in WebAssembly)

> [!IMPORTANT]
> **Investigation required.** Measure web2 with the current build against the CLI on `data/zika/20` and on the datasets with the largest earlier changes (`ebola/20` with coalescent optimization, `flu/h3n2/20` dense, `lassa/L/50`), and find the steps of `timetree` that amplify rounding differences when the difference exceeds the threshold. The core test suite can run in WebAssembly under Node (`cargo test --target wasm32-unknown-emscripten` with unwinding and a Node runner), which compares each step with its native result within the test tolerances.

The reproducibility target is equality within tolerances ([kb/decisions/reproducibility-across-platforms-within-tolerances.md](../decisions/reproducibility-across-platforms-within-tolerances.md)): the steps of `timetree` that amplify rounding differences beyond the threshold must be made stable.

## Validation

- The clock rate, the node dates, and the trees of `timetree` agree between the CLI with different OpenBLAS kernels and web2 within the agreed tolerance, or byte for byte
