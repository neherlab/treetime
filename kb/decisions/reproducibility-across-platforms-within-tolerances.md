# Outputs agree across platforms within tolerances

TreeTime run on different CPUs, with different OpenBLAS builds, natively and in the WebAssembly back end of web2, gives outputs that agree within a tolerance per output. Byte-identical outputs across platforms are not a goal.

**Type**: Design rule for v1.

**Status**: Approved on 2026-10-21.

## Rule

- **Continuous outputs** (rates, branch lengths, dates, probabilities, likelihoods) agree within tolerances per output. A step that amplifies rounding differences beyond the project threshold of 1e-6 relative is a defect to fix in the algorithm
- **Discrete outputs** taken from continuous values (a reconstructed state or character, a date on a grid, a root edge) follow a rule that gives the same choice on every platform when the values are equal up to rounding
- **Platforms keep their fastest arithmetic**: OpenBLAS selects its kernels per CPU (with or without FMA), and the math library and LAPACK build may differ between targets

## Why

Byte-identical outputs need one BLAS kernel set without FMA, one LAPACK build, one math library, and fixed kernels for small products and FFT on every target. That costs the speed of the FMA kernels on every CPU, and a build against a system OpenBLAS loses the guarantee anyway. Agreement within tolerances keeps the fast kernels and asks the algorithms to be numerically stable, which their results need regardless of the platform.

## Related

- [kb/issues/M-timetree-outputs-depend-on-the-rounding-of-the-platform.md](../issues/M-timetree-outputs-depend-on-the-rounding-of-the-platform.md): measured differences between platforms
- [kb/issues/M-posterior-argmax-follows-rounding-at-ties.md](../issues/M-posterior-argmax-follows-rounding-at-ties.md): a tied state follows the rounding of the platform
