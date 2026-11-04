# The FFT round-off bound assumes that `sin` and `cos` are within 1 ulp

> [!IMPORTANT]
> **Investigation required.** No math library that TreeTime runs on documents an error bound for `sin` and `cos`, so the premise of the FFT round-off bound is unproven. Evidence to collect: a proof for every platform, which correctly rounded `sin` and `cos` in every artifact give by construction

## Problem

`fn FftConvolver::overlap_sums_round_off()` ([packages/treetime-ops/src/convolution.rs#L173-L194](../../packages/treetime-ops/src/convolution.rs#L173-L194)) bounds the round-off of an FFT convolution by the error analysis of [kb/decisions/distribution-convolution-tails-by-exponential-tilting.md](../decisions/distribution-convolution-tails-by-exponential-tilting.md) (section "Exact tails at hard ends"). Its twiddle term $\mu = 2\pi\gamma_2 + \sqrt{2}\,u$ assumes that each computed twiddle component is within 1 ulp of the cosine or sine of the computed angle. The windowed tilts of the exact tails keep a sample only where this bound is at most $10^{-6}$ nats, and the exact-tail tests assert the bound, so the premise decides which samples the FFT computes.

`rustfft` 6.4.1 (`fn compute_twiddle()` in `src/twiddles.rs`) and `realfft` 3.5.0 (`fn compute_twiddle()` in `src/lib.rs`) compute every twiddle as `f64::cos` and `f64::sin` of the angle $(-2\pi/N) \cdot k$, with $|\theta| < 2\pi$. The Rust standard library calls the math library of the platform: glibc or musl on Linux, the Apple library on macOS, the Microsoft C runtime on Windows, and the musl-based library of Emscripten in the WebAssembly builds.

## Evidence

The largest errors known for binary64 `sin` and `cos` (the larger of the two), found by search over all inputs (Gladman, Innocente, Mather, Ozaki, Zimmermann, *Accuracy of Mathematical Functions in Single, Double, Double Extended, and Quadruple Precision*, August 2026, Table 3, https://members.loria.fr/PZimmermann/papers/accuracy.pdf):

- glibc 2.44: 0.516 ulp
- musl 1.2.6: 0.834 ulp
- Apple 26.5.2: 0.948 ulp
- Microsoft (MSVC 2026): 0.897 ulp
- LLVM libc 22.1.8: 0.500 ulp (correctly rounded)

The paper states that these values are lower bounds of the true largest errors. Every known value is below 1 ulp, but no library proves it.

## Sensitivity

With errors up to $k$ ulp, the twiddle term becomes $\mu = 2\pi\gamma_2 + \sqrt{2}\,k u$. Because $2\pi\gamma_2 \approx 12.6u$ dominates $\mu$, the level error $\eta$ grows by 7% for $k = 2$ and by 21% for $k = 4$, and the bound $\delta$ by about the same factor. A bound that is too small by such a factor keeps samples whose round-off reaches about $1.07 \cdot 10^{-6}$ or $1.2 \cdot 10^{-6}$ nats instead of $10^{-6}$.

## Fix direction

Link correctly rounded `sin` and `cos` (for example those of LLVM libc, within 0.5 ulp by construction) into every artifact, so that the premise holds on every platform. Then state the premise as proven in the decision file and delete this issue. A test that compares every twiddle with a multiple-precision reference (MPFR) up to the largest FFT length proves the premise only for the math library of the host that runs the test.
