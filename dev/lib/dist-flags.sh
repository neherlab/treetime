#!/usr/bin/env bash
#
# The CPU flags and the unwinder link of the shipped builds, shared by the cross
# builds (dev/cross/build) and the native dist, profiling, bench, and addon
# builds of the justfile, so a local `just build-dist` compiles and links the
# code that ships. For x86_64-unknown-linux-gnu these are all the code
# generation and link flags of the shipped binary; dev/cross/build adds the
# static runtime and linker flags of the other targets.
#
# The prebuilt OpenBLAS, liblzma, and libzstd archives that the builds link are
# compiled with these C CPU flags; the images check the build record of every
# archive against them (dev/docker/files/install-image-libraries). OpenBLAS also
# selects its kernels at run time.

# Usage: export_dist_flags <binary|addon> [target]
#
# Export RUSTFLAGS and the C and C++ flags of the shipped build of a binary or
# of the Node addon for a target, by default the host. RUSTFLAGS replaces the
# target rustflags of every Cargo config file (frame pointers and v0 symbol
# mangling of .cargo/config.toml, mold of the development image) instead of
# adding to them, so the native builds link with rust-lld. The C and C++ flags
# carry only the CPU baseline: the cc crate passes the optimization level of
# the cargo profile.
export_dist_flags() {
  local kind="${1:?binary or addon}" target="${2:-$(rustc --print host-tuple)}" target_safe c_arch rust_arch link=""
  target_safe="${target//-/_}"

  case "${kind}" in
  binary | addon) ;;
  *)
    printf 'Unknown build kind: %s. Supported: binary, addon\n' "${kind}" >&2
    return 1
    ;;
  esac

  # The glibc binaries link the GCC unwinder statically, so they need no
  # libgcc_s.so.1 and run on images that carry glibc alone (see the linker
  # script in dev/cross/static-unwinder/). The Node addon keeps the shared
  # unwinder of the Electron process that loads it.
  if [[ "${kind}" == binary && "${target}" == *-linux-gnu ]]; then
    link=" -L native=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)/dev/cross/static-unwinder"
  fi

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

  export RUSTFLAGS="${rust_arch}${link}"
  export "CFLAGS_${target_safe}=${CFLAGS:-} ${c_arch}"
  export "CXXFLAGS_${target_safe}=${CXXFLAGS:-} ${c_arch}"
}
