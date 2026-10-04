#!/usr/bin/env bash
#
# The CPU flags of the shipped builds, shared by the cross builds
# (dev/cross/build) and the native dist, profiling, and bench builds of the
# justfile, so a local `just build-dist` compiles the code that ships. For
# x86_64-unknown-linux-gnu these are all the code generation flags of the
# shipped binary; dev/cross/build adds the static runtime and linker flags of
# the other targets.

# Export RUSTFLAGS and the C, C++, and Fortran flags of the shipped build for a
# target, by default the host. RUSTFLAGS replaces the target rustflags of
# .cargo/config.toml (mold, frame pointers, v0 symbol mangling, OpenBLAS link
# arguments) instead of adding to them, as in the shipped builds.
export_dist_flags() {
  local target="${1:-$(rustc --print host-tuple)}" target_safe c_arch rust_arch c_flags
  target_safe="${target//-/_}"

  case "${target}" in
  aarch64-apple-*)
    c_arch="-mcpu=apple-m1"
    rust_arch="-C target-cpu=apple-m1"
    ;;
  aarch64-*)
    # Features of -C target-feature=+v8.2a, which is unstable, listed explicitly.
    c_arch="-march=armv8.2-a"
    rust_arch="-C target-cpu=generic -C target-feature=+crc,+dpb,+lor,+lse,+neon,+pan,+ras,+rdm,+vh"
    ;;
  x86_64-*)
    c_arch="-march=haswell"
    rust_arch="-C target-cpu=haswell"
    ;;
  *)
    printf 'Unsupported target: %s\n' "${target}" >&2
    return 1
    ;;
  esac

  c_flags="-O3 -ftree-vectorize -funroll-loops ${c_arch}"
  export RUSTFLAGS="${rust_arch}"
  export "CFLAGS_${target_safe}=${CFLAGS:-} ${c_flags}"
  export "CXXFLAGS_${target_safe}=${CXXFLAGS:-} ${c_flags}"
  export "FCFLAGS_${target_safe}=${FCFLAGS:-} ${c_flags}"
}
