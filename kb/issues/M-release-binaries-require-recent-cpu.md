# Release binaries require a recent CPU

## Summary

`export_dist_flags()` in [dev/lib/dist-flags.sh](../../dev/lib/dist-flags.sh) sets the CPU baseline of every shipped build (`dev/cross/build`, `just build-dist`, `just build-napi dist`): Haswell (x86-64-v3 with AVX2 and FMA) on x86_64, and ARMv8.2-A extensions (LSE atomics, RDM, and others) on Linux and Windows aarch64. The C dependencies get the same baseline through `-march`. The binaries do not run on an older CPU. The developer guide states the requirement, but no user-facing document, release note, or CLI help does.

## Evidence

- `dev/lib/dist-flags.sh` sets `-C target-cpu=haswell` and `-march=haswell` for every `x86_64*` target except `x86_64-win7-*`, which gets `x86-64-v2`. The Windows 7 target is not in `dev/cross/targets`, so the nightly releases do not ship it.
- `dev/lib/dist-flags.sh` sets `-C target-cpu=generic -C target-feature=+crc,+dpb,+lor,+lse,+neon,+pan,+ras,+rdm,+vh` and `-march=armv8.2-a` for non-Apple `aarch64*` targets.
- On a CPU without the build features, the CLI and the server stop at startup with a message that names the required CPU level: `require_cpu_features()` [packages/treetime-utils/src/init/cpu.rs#L5-L14](../../packages/treetime-utils/src/init/cpu.rs#L5-L14), called first in [packages/app-cli/src/bin/treetime.rs#L18](../../packages/app-cli/src/bin/treetime.rs#L18) and [packages/app-server/src/bin/treetime-server.rs#L34](../../packages/app-server/src/bin/treetime-server.rs#L34). The Node addon returns the same message as an error ([packages/app-napi/src/guard.rs#L8-L9](../../packages/app-napi/src/guard.rs#L8-L9)).
- [docs/dev/developer_guide.md#L443](../../docs/dev/developer_guide.md#L443), section "Cross-compilation and releases", states the requirement for developers.
- TreeTime v0 is a Python package and runs on any CPU that runs NumPy.

## Impact

- The binaries do not run on x86_64 CPUs older than Haswell (2013), which are still common in HPC clusters and older servers.
- The Linux and Windows aarch64 binaries do not run on ARMv8.0 cores such as the Cortex-A53 and Cortex-A72 (Raspberry Pi 3 and 4, AWS Graviton1).
- Users learn the requirement only when the binary refuses to start.

## Potential solutions

> [!IMPORTANT]
> **Decision required.** The CPU baseline of the shipped binaries is open. Options:
>
> - Lower the baseline to `x86-64-v2` and generic `aarch64`. This costs some performance in vectorized code and in OpenBLAS.
> - Ship two variants per architecture (baseline and optimized) and select one at install time or with a small launcher. This doubles the release artifacts.
> - Keep the baseline and state the CPU requirement in the release notes and the installation instructions.
>
> Evidence: the startup check already turns the failure into a clear message, so the remaining cost of the current baseline is the excluded hardware, not a silent crash. No measurement of the performance difference between the baselines exists.
