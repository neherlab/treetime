# The outputs of timetree depend on the rounding of the platform

Differences in the last bits of linear algebra and math functions grow through the time inference of `timetree` to relative differences of about 2e-5 in the clock rate, above the project threshold of 1e-6. The same inputs and configuration give different clock rates, node dates, and trees on the CLI with the CPU-specific kernels of OpenBLAS, on the CLI with its generic kernels, and in the WebAssembly back end of web2. `ancestral` gives byte-identical outputs on these platforms.

## Evidence

Measured on `data/zika/20` (`timetree` with the configuration that web2 resolves for a run with tree, alignment, and metadata: `max_iter` 2, marginal branch lengths, `clock_filter` 3), single-threaded, each platform reproducible run to run:

> [!WARNING]
> **Needs review.** These values were measured with OpenBLAS 0.3.28 natively and OpenBLAS develop `eb00d597` in web2. The builds now use OpenBLAS develop `de5d4914` on every target: natively with LAPACK compiled by flang 22, in web2 with the LAPACK C translation (`dev/docker/files/install-openblas`, `dev/docker/files/install-openblas-emscripten`), so the native values differ; the LAPACK C translation of web2 is the same in both develop commits

| Output                                                  | CLI, OpenBLAS `Haswell` kernels | CLI, `OPENBLAS_CORETYPE=Prescott` (`Katmai` kernels) | web2 (Emscripten, `WASM128_GENERIC`) |
| ------------------------------------------------------- | ------------------------------- | ---------------------------------------------------- | ------------------------------------ |
| `clock_rate` of `timetree.clock-model.json`             | 0.0010266268996163546           | 0.0010266466684456434                                | 0.001026646778721711                 |
| `log_lh_seq` of the first iteration (`tracelog.csv`)    | -17281.72723068383              | -17281.72723068391                                   | -17281.73070320649                   |
| `log_lh_pos` of the first iteration                     | -92.15030048682337              | -91.88381677167982                                   | -92.14962165070715                   |
| `timetree.gtr.json`, `timetree.reconstructed-nuc.fasta` | reference                       | identical                                            | identical                            |

- The clock rate of web2 differs from the CLI with `Haswell` kernels by 1.9e-5 relative, and from the CLI with `Katmai` kernels by 1.1e-7 relative: the spread between two CPU kernel sets of the native build is as large as the difference between web2 and the CLI
- The root-to-tip divergence of the first tip differs between web2 and the CLI by 1.2e-9 relative; the log-likelihoods of the first iteration differ by 2.0e-7 relative (sequences) and by 0.27 absolute between the two native kernel sets (positions), so the time inference amplifies differences far below the threshold
- The OpenBLAS version alone changes `timetree` beyond the threshold. `dev/smoke` with the same code, once with OpenBLAS 0.3.28 (commit `2dabfc59`) and once with 0.3.34, both with `Haswell` kernels on an Intel i9-10900K: the clock rate changes by 9.2e-4 relative on `lassa/L/50`, 3.9e-4 on `ebola/362` with skyline, 3.2e-4 on `ebola/100`, 2.6e-4 on `tb/20`, and 7.9e-6 on `rsv/a/100`; the branch lengths of the time tree change by up to 0.4 on `lassa/L/50`, and the reconstructed sequences of `lassa/L/50` and of `ebola/362` with `--clock-filter=10` change. `timetree` on `zika/20` stays byte-identical. The outputs of `clock`, `homoplasy`, and `prune` stay byte-identical, and those of `ancestral`, `mugration`, and `optimize` change by at most 1.8e-7 relative, apart from the last printed digit of one value in `mugration.confidence.csv`
- The OpenBLAS build changes `timetree` beyond the threshold again. `dev/smoke` (quick tier) of commit `2e275920` against `4f116c6d`, the same code once with the OpenBLAS 0.3.34 release of the earlier images, once with develop `de5d4914` and LAPACK compiled by flang 22, both with `Haswell` kernels on an Intel i9-10900K: `timetree` changes in 58 cases; the clock rate changes by 1.7e-3 relative on `ebola/20` with coalescent optimization, 1.5e-3 on `flu/h3n2/20` dense and with coalescent, 3.7e-4 on `dengue/20`, 3.4e-4 on `ebola/362` with skyline, 3.3e-4 on `ebola/100`, 3.2e-4 on `tb/20`, below 3e-4 on the other datasets, and by at most 2.4e-11 on most `flu/h3n2/20` variants. The other commands change by at most: `ancestral` 1.9e-12 relative (GTR of `--gtr-iterations=3`), `optimize` 4.3e-11 apart from one branch of 3.7e-8 on `ebola/100` (1.1e-5 relative, 4.2e-13 absolute), `mugration` 1.4e-7 apart from the last printed digit of one value in `mugration.confidence.csv` and one tied state on `dengue/100` ([N-mugration-tied-state-follows-rounding.md](N-mugration-tied-state-follows-rounding.md)). `clock`, `homoplasy`, and `prune` stay byte-identical. `cli/pipeline-steps-timetree` and `config/zika/20/pipeline` change through their `timetree` step
- `ancestral` on `data/flu/h3n2/20`: `augur-node-data.json`, `gtr.json`, `nexus`, `nwk`, and `reconstructed-nuc.fasta` of web2 and of the CLI are byte-identical
- Rounding sources that differ between the platforms: the OpenBLAS kernels (FMA in `Haswell`, none in `Katmai` and `WASM128_GENERIC`); LAPACK (OpenBLAS develop `de5d4914` compiled by flang 22 natively, its C translation with the LAPACK 3.9 `?lartg` and `?lassq` in web2); `exp`, `log`, and `pow` (glibc natively, with FMA variants on FMA CPUs; musl of Emscripten, from the same ARM optimized routines, in web2); the small products of `matrixmultiply` (FMA kernels natively, the fallback kernel in WebAssembly)

> [!IMPORTANT]
> **Investigation required.** Find the steps of `timetree` that amplify rounding differences beyond the threshold (the positional likelihood of the first iteration differs by 0.27 between two native kernel sets), and which of the sources above causes the difference in the sequence likelihood between web2 and the CLI. The core test suite can run in WebAssembly under Node (`cargo test --target wasm32-unknown-emscripten` with unwinding and a Node runner), which compares each step with its native result within the test tolerances.

> [!IMPORTANT]
> **Decision required.** The reproducibility target across platforms:
>
> - equality within tolerances: the CLI on every CPU and web2 agree within agreed tolerances per output, and the time inference is made stable enough that rounding differences stay within them
> - byte-identical outputs: one BLAS kernel set without FMA, one LAPACK build, one math library, and fixed kernels for small products and FFT on every target, natively and in web2

## Validation

- The clock rate, the node dates, and the trees of `timetree` agree between the CLI with different OpenBLAS kernels and web2 within the agreed tolerance, or byte for byte
