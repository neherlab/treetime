# The desktop build drops the C include path the native addon needs

## Summary

`just build-desktop` builds the N-API addon first (`@neherlab/app-napi#build`, through turbo). Turbo runs tasks in strict environment mode and passes only the variables listed in `passThroughEnv`. The build container sets `C_INCLUDE_PATH=/usr/include/x86_64-linux-gnu` and `CPLUS_INCLUDE_PATH` so that its GCC finds the Debian multiarch headers, but the addon task does not pass them through. The C build scripts of `zstd-sys` and `lzma-sys` then fail, and the desktop build stops before the Electron bundles are built.

## Evidence

- The addon task passes `CARGO_*`, `RUSTUP_*`, `KACHE_*`, `PATH`, `HOME`, `LIBRARY_PATH`, `LD_LIBRARY_PATH` and `PKG_CONFIG_PATH`, but not `C_INCLUDE_PATH` or `CPLUS_INCLUDE_PATH` [`packages/app-napi/turbo.json`](../../packages/app-napi/turbo.json)
- In the container, `just build-desktop` fails with `/usr/include/stdlib.h:26:10: fatal error: bits/libc-header-start.h: No such file or directory` while compiling `zstd/lib/legacy/zstd_v07.c` and `xz-5.2/src/liblzma/simple/x86.c`
- The same `cc` compiles `zstd_v07.c` without error when run directly in the container, where `C_INCLUDE_PATH` is set

## Impact

- `just build-desktop` cannot produce the desktop app in the build container
- The renderer, main and preload bundles build with `vite build` in `packages/app-desktop`, so front-end changes can still be checked
