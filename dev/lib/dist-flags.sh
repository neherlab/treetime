#!/usr/bin/env bash
#
# The CPU flags of the shipped builds, shared by the cross builds
# (dev/cross/build) and the native dist, profiling, and bench builds of the
# justfile, so a local `just build-dist` compiles the code that ships. For
# x86_64-unknown-linux-gnu these are all the code generation flags of the
# shipped binary; dev/cross/build adds the static runtime and linker flags of
# the other targets.
#
# The prebuilt archives that the builds link follow the same CPU baselines:
# OpenBLAS selects its kernels at run time, and binarylandia/staticbourg builds
# liblzma and libzstd for these baselines (.config/mise.image.toml).

# Export RUSTFLAGS and the C and C++ flags of the shipped build for a target, by
# default the host. RUSTFLAGS replaces the target rustflags of .cargo/config.toml
# (mold, frame pointers, v0 symbol mangling) instead of adding to them, as in
# the shipped builds. The C and C++ flags carry only the CPU baseline: the cc
# crate passes the optimization level of the cargo profile.
export_dist_flags() {
  local target="${1:-$(rustc --print host-tuple)}" target_safe c_arch rust_arch
  target_safe="${target//-/_}"

  case "${target}" in
  aarch64-apple-*)
    c_arch="-mcpu=apple-m1"
    rust_arch="-C target-cpu=apple-m1"
    ;;
  aarch64-*-linux-*)
    # Armv8.0, so the binaries also run on Cortex-A53 and Cortex-A72 cores
    # (Raspberry Pi 3 and 4). Atomic operations go through the outline atomics
    # helpers, which use LSE instructions when the CPU has them; rustc enables
    # them by default for these targets, and -moutline-atomics does the same for
    # the C dependencies.
    c_arch="-march=armv8-a -moutline-atomics"
    rust_arch="-C target-cpu=generic"
    # jemalloc fixes its page size at build time and aborts at startup on a
    # kernel with larger pages. 64K pages also work on kernels with 4K and 16K
    # pages (Raspberry Pi 5, Asahi Linux, Arm servers with 64K kernels).
    export JEMALLOC_SYS_WITH_LG_PAGE=16
    ;;
  aarch64-*)
    # Windows 11 on Arm requires the LSE atomics of Armv8.1, and the Windows
    # targets of rustc have no outline atomics, so LSE is enabled at build time.
    c_arch="-march=armv8-a+lse"
    rust_arch="-C target-feature=+lse"
    ;;
  x86_64-*)
    # x86-64-v2 (SSE4.2, POPCNT, CMPXCHG16B; Intel Nehalem, AMD Bulldozer, or
    # newer), so the binaries also run on CPUs and virtual machines without AVX2.
    c_arch="-march=x86-64-v2"
    rust_arch="-C target-cpu=x86-64-v2"
    ;;
  *)
    printf 'Unsupported target: %s\n' "${target}" >&2
    return 1
    ;;
  esac

  export RUSTFLAGS="${rust_arch}"
  export "CFLAGS_${target_safe}=${CFLAGS:-} ${c_arch}"
  export "CXXFLAGS_${target_safe}=${CXXFLAGS:-} ${c_arch}"
}
