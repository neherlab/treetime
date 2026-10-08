# Web2 is not verified in browsers

The WebAssembly back end of web2 has run only in Node and Bun (`packages/app-wasm/smoke.ts`, and probes that send API requests to `Backend.fetch()`). The page, the worker start, the message port transport, and the downloads have not run in Chrome, Firefox, or Safari.

## Evidence

- `just test-wasm` runs the smoke test in Bun; TypeScript tests of the project run without browsers by project rule
- Browser behavior that Node does not show: the content security policy of the page (`'wasm-unsafe-eval'`), the memory limit per page (iOS Safari allows much less than 4 GiB, [emscripten#19144](https://github.com/emscripten-core/emscripten/issues/19144)), and transfer of large response chunks over the message port

> [!IMPORTANT]
> **Investigation required.** Run web2 in current Chrome, Firefox, and Safari: start the back end, upload a dataset, run `ancestral` and `timetree`, download a file and the run archive, and record the worker memory.

## Validation

- The steps above pass in each browser that web2 supports, and the supported browsers are documented
