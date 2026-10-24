# Web2 is not verified in Firefox and Safari

The WebAssembly back end of web2 has run in Node, Bun, and Chrome. The page, the worker start, the message port transport, the threaded module, and the downloads have not run in Firefox or Safari, and the memory of the workers has not been recorded in any browser.

## Evidence

- `just test-wasm` and `just test-wasm-mt` run the smoke test in Bun and Node; TypeScript tests of the project run without browsers by project rule
- Chrome on Linux runs the production build and the dev server of web2: the page is cross-origin isolated, an example dataset loads from the examples panel and from "Example files", a `sc2/5100` timetree run shows live progress and finishes in about 11 s on the threaded module and 23 s with `?jobs=1`, Stop ends a run as cancelled, a run file and the run archive download, and a page over plain HTTP on the LAN, which is not isolated, finishes runs on the single-threaded module
- Browser behavior that Node does not show: the content security policy of the page (`'wasm-unsafe-eval'`), cross-origin isolation and shared memory, the memory limit per page (iOS Safari allows much less than 4 GiB, [emscripten#19144](https://github.com/emscripten-core/emscripten/issues/19144)), and transfer of large response chunks over the message port

> [!IMPORTANT]
> **Investigation required.** Run web2 in current Firefox and Safari: start the back end, load an example dataset, run `ancestral` and `timetree` with and without `?jobs=1`, stop a run, download a file and the run archive, and record the memory of the back-end and job workers, also in Chrome.

## Validation

- The steps above pass in each browser that web2 supports, and the supported browsers are documented
