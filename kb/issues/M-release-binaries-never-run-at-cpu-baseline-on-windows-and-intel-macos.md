# Binaries never run at the CPU baseline on Windows and Intel macOS

## Summary

The Windows x86_64 binaries (UCRT and Windows 7), the Windows arm64 binaries, and the macOS x86_64 binaries never run on a CPU at their baseline. CI runs them only on runners above it, or not at all (Windows 7), so an unguarded instruction above the baseline in TreeTime or in a dependency shows only on a user's machine.

## Evidence

- Baselines (`export_dist_flags` in `dev/lib/dist-flags.sh`): x86-64-v2 on every x86_64 target, Armv8.0 with LSE on Windows arm64, Armv8.0 on Linux arm64, and the Apple M1 on macOS arm64
- Linux binaries run under QEMU on the oldest CPU of their baseline (Nehalem, Cortex-A53), and on x86_64 also on a CPU below it (`dev/cross/check-cpus`, CI job `cpu-baseline` in `.github/workflows/cli.yml`, and `dev/cross/all --run`)
- macOS arm64 binaries run on Apple Silicon runners (`compat-macos`), whose oldest model is the baseline, the M1
- Windows binaries run on `windows-2022` and `windows-11-arm` (`compat-windows`), and macOS x86_64 binaries on `macos-15-intel` and `macos-26-intel` (`compat-macos`). These runners have CPUs far above x86-64-v2 and Armv8.0 with LSE
- The Windows 7 binary (`x86_64-win7-windows-gnu`) is not in `dev/cross/targets` and runs in no CI job
- The Node addons of the desktop app (`app-napi`) run on a baseline CPU on no platform; they link the same code and prebuilt archives as the CLI

Existing safeguards:

- `export_dist_flags` owns the CPU flags of every shipped build, Rust and C
- The image builds check the build record of every prebuilt archive (OpenBLAS, liblzma, libzstd) against these flags (`dev/docker/files/install-image-libraries`)

Neither safeguard sees code that a dependency compiles with its own flags for a newer CPU and calls without a check of the CPU at run time.

> [!IMPORTANT]
> **Investigation required.** Check the license terms of Intel SDE for use in CI, a pinned download with a checksum, whether the current 10.x series still runs on macOS, and the run time of the `dev/cross/check-cpus` workload under SDE.

## Impact

- A release can crash with an illegal instruction on an older Windows or Intel Mac that the documented baseline supports, and CI passes

## Potential solutions

- Run the x86_64 Windows binaries and the Windows 7 binary on a Windows runner under Intel SDE with its chip check and the CPUID of the baseline CPU (`sde64 -nhm`, Nehalem for x86-64-v2), with the workload of `dev/cross/check-cpus`
- Run the macOS x86_64 binary under SDE in the same way, if SDE still runs on macOS
- Windows arm64 has no user-mode emulator for an older CPU, so its option stays open

## Sources

- https://www.intel.com/content/www/us/en/developer/articles/tool/software-development-emulator.html
- https://github.com/jvdd/argminmax/pull/87 (SDE in a Rust CI job)
