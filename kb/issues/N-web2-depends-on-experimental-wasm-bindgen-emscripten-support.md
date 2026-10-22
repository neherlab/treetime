# Web2 depends on the experimental Emscripten support of wasm-bindgen

`emcc` runs the wasm-bindgen CLI on the linked module (`-sWASM_BINDGEN` in `[target.wasm32-unknown-emscripten]` of `.cargo/config.toml`). Both projects mark this integration experimental, so an upgrade of either can change flags or output and break the web2 build or the generated declarations.

## Evidence

- The wasm-bindgen guide (`guide/src/reference/emscripten.md`) says that Emscripten support "is still being smoothed out; flags and output shape may still change". It came with wasm-bindgen 0.2.115 (2026-03-27, [wasm-bindgen#4443](https://github.com/wasm-bindgen/wasm-bindgen/pull/4443))
- Emscripten marks the `WASM_BINDGEN` setting experimental since 6.0.9; the post-link flow that web2 uses needs 6.0.10 (Emscripten `ChangeLog.md`)
- Unwinding is not supported on this target ([wasm-bindgen#5165](https://github.com/wasm-bindgen/wasm-bindgen/issues/5165)), see [M-web2-wasm-panic-loses-the-back-end-state.md](M-web2-wasm-panic-loses-the-back-end-state.md)
- The threaded module of the job workers depends on the patched build of the wasm-bindgen CLI from `binarylandia/appville` (task `tasks/wasm-bindgen`), which both builds use: upstream 0.2.129 moves the start function of the linker, which initializes the shared memory of a threaded module, to `__wbindgen_start`, which its JavaScript calls after instantiation. The memory is then no longer initialized at instantiation, and the `emcc` link of a threaded module whose `EM_JS` or `EM_ASM` strings `emcc` reads from that memory fails (for example with `-sMALLOC=mimalloc`). An upgrade of wasm-bindgen therefore also needs the patch to apply to the new version and a new appville release
- The interface of the module that wasm-bindgen generates is small: `Backend`, `Exchange`, `ResponseHead`, `ChunkReader`, `JobOrder`, and `Job` in `packages/app-wasm/src/`

> [!IMPORTANT]
> **Decision required.** How web2 handles the experimental status:
>
> - keep wasm-bindgen and upgrade wasm-bindgen and Emscripten together, checked by a CI job ([N-web2-production-build-not-in-ci.md](N-web2-production-build-not-in-ci.md))
> - replace wasm-bindgen with C exports of the module that the worker calls through the Emscripten runtime, with hand-written declarations for the few exported functions

## Validation

- The chosen approach builds and passes `just test-wasm` with the pinned toolchain
