# Release binaries require a recent CPU

## Summary

`dev/cross/build` compiles the release binaries for a fixed CPU baseline: Haswell (x86-64-v3 with AVX2 and FMA) on x86_64, and ARMv8.2-A extensions (LSE atomics, RDM, and others) on Linux aarch64. The C dependencies get the same baseline through `-march`. A binary built this way stops with `SIGILL` (illegal instruction) on an older CPU. Nothing in the release notes, the CLI help, or the documentation states the requirement.

## Evidence

- `dev/cross/build` sets `-C target-cpu=haswell` and `-march=haswell` for every `x86_64*` target.
- `dev/cross/build` sets `-C target-cpu=generic -C target-feature=+crc,+dpb,+lor,+lse,+neon,+pan,+ras,+rdm,+vh` and `-march=armv8.2-a` for `aarch64*` Linux targets.
- TreeTime v0 is a Python package and runs on any CPU that runs NumPy.

## Impact

- The binaries fail on x86_64 CPUs older than Haswell (2013), which are still common in HPC clusters and older servers.
- The Linux aarch64 binaries fail on ARMv8.0 cores such as the Cortex-A53 and Cortex-A72 (Raspberry Pi 3 and 4, AWS Graviton1).
- The failure is a crash with no error message, so users cannot tell the cause.

## Potential solutions

- Lower the baseline to `x86-64-v2` and generic `aarch64`. This costs some performance in vectorized code.
- Ship two variants per architecture (baseline and optimized) and select one at install time or with a small launcher.
- Keep the baseline and document the CPU requirement in the release notes and the installation instructions, and print a clear error at startup when the CPU lacks a required feature.
