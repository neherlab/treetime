# Release binaries require a recent CPU

## Summary

`export_dist_flags()` in [dev/lib/dist-flags.sh](../../dev/lib/dist-flags.sh) sets the CPU baseline of every shipped build (`dev/cross/build`, `just build-dist`, `just build-napi dist`): Haswell (x86-64-v3 with AVX2 and FMA) on x86_64, and ARMv8.2-A extensions (LSE atomics, RDM, and others) on Windows aarch64. The C dependencies get the same baseline through `-march`. The x86_64 binaries do not run on an older CPU. The developer guide states the requirement, but no user-facing document, release note, or CLI help does.

The Linux aarch64 binaries are built for ARMv8.0 with outline atomics and run on every 64-bit Arm CPU.

## Evidence

- `dev/lib/dist-flags.sh` sets `-C target-cpu=haswell` and `-march=haswell` for every `x86_64*` target except `x86_64-win7-*`, which gets `x86-64-v2`. The Windows 7 target is not in `dev/cross/targets`, so the nightly releases do not ship it.
- `dev/lib/dist-flags.sh` sets `-C target-cpu=generic -C target-feature=+crc,+dpb,+lor,+lse,+neon,+pan,+ras,+rdm,+vh` and `-march=armv8.2-a` for `aarch64-pc-windows-gnullvm`. Windows 11 on Arm requires ARMv8.1 atomics, and every Windows 11 on Arm device has an ARMv8.2 core, so this baseline excludes only Windows 10 devices with a Snapdragon 835.
- On a CPU without the build features, the x86_64 CLI and server stop at startup with a message that names the required CPU level: `require_cpu_features()` [packages/treetime-utils/src/init/cpu.rs#L5-L14](../../packages/treetime-utils/src/init/cpu.rs#L5-L14), called first in [packages/app-cli/src/bin/treetime.rs#L18](../../packages/app-cli/src/bin/treetime.rs#L18) and [packages/app-server/src/bin/treetime-server.rs#L34](../../packages/app-server/src/bin/treetime-server.rs#L34). The Node addon returns the same message as an error ([packages/app-napi/src/guard.rs#L8-L9](../../packages/app-napi/src/guard.rs#L8-L9)). The Windows aarch64 build has no such check.
- [docs/dev/developer_guide.md#L443](../../docs/dev/developer_guide.md#L443), section "Cross-compilation and releases", states the requirements for developers.
- TreeTime v0 is a Python package and runs on any CPU that runs NumPy.
- A lower baseline cost no measurable time on aarch64: on a GitHub `ubuntu-24.04-arm` runner (4 Neoverse N2 cores), the ARMv8.0 build with outline atomics and the ARMv8.2 build ran the `flu/h3n2/500` pipeline (4 threads and 1 thread) and the dense marginal reconstruction of `dengue/2000` within 0.5% of each other. The same runs with `OPENBLAS_CORETYPE` set to `ARMV8`, `NEOVERSEN1`, and `NEOVERSEN2` were also within 0.5%. No such measurement exists for x86_64, where the baseline decides whether the compiler can use AVX2 and FMA.

## Impact

- The x86_64 binaries do not run on CPUs older than Haswell (2013), which are still common in HPC clusters and older servers.
- Users learn the requirement only when the binary refuses to start.

## Potential solutions

> [!IMPORTANT]
> **Decision required.** The x86_64 CPU baseline of the shipped binaries is open. Options:
>
> - Lower the baseline to `x86-64-v2`. This can cost performance in code that the compiler vectorizes with AVX2 and FMA; OpenBLAS selects its kernels at run time and is not affected.
> - Ship two variants (baseline and optimized) and select one at install time or with a small launcher. This doubles the x86_64 release artifacts.
> - Keep the baseline and state the CPU requirement in the release notes and the installation instructions.
>
> Evidence: the startup check already turns the failure into a clear message, so the remaining cost of the current baseline is the excluded hardware, not a silent crash. The aarch64 measurement above does not transfer to x86_64; an x86_64 measurement of `x86-64-v2` against `haswell` with the same commands is missing.
