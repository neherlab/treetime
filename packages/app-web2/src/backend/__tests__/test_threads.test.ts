import { chooseRunThreads, DEFAULT_JOBS, parseJobs, type ThreadScope } from "../threads";

const THREADED_WASM_URL = "https://treetime.localhost/assets/app_wasm-mt.wasm";

describe("jobs URL parameter", () => {
  test.each([
    ["no parameter", "", { jobs: DEFAULT_JOBS, warning: undefined }],
    ["one thread", "?jobs=1", { jobs: 1, warning: undefined }],
    ["four threads", "?example=flu&jobs=4", { jobs: 4, warning: undefined }],
    [
      "zero",
      "?jobs=0",
      { jobs: DEFAULT_JOBS, warning: "the URL parameter jobs=0 is not a positive integer; runs use 8 threads" },
    ],
    [
      "a word",
      "?jobs=many",
      { jobs: DEFAULT_JOBS, warning: "the URL parameter jobs=many is not a positive integer; runs use 8 threads" },
    ],
    [
      "a fraction",
      "?jobs=2.5",
      { jobs: DEFAULT_JOBS, warning: "the URL parameter jobs=2.5 is not a positive integer; runs use 8 threads" },
    ],
  ])("%s", (_name, search, expected) => {
    expect(parseJobs(search)).toStrictEqual(expected);
  });
});

describe("run module choice", () => {
  test("one thread selects the single-threaded module without probing the browser", () => {
    expect(chooseRunThreads(1, scope({ crossOriginIsolated: false }), THREADED_WASM_URL)).toStrictEqual({
      jobs: 1,
      threadedWasmUrl: undefined,
      warning: undefined,
    });
  });

  test("several threads in a cross-origin isolated page with shared memory select the threaded module", () => {
    expect(chooseRunThreads(4, scope(), THREADED_WASM_URL)).toStrictEqual({
      jobs: 4,
      threadedWasmUrl: THREADED_WASM_URL,
      warning: undefined,
    });
  });

  test.each([
    [
      "a page that is not cross-origin isolated",
      { crossOriginIsolated: false },
      "the page is not cross-origin isolated",
    ],
    ["a browser without SharedArrayBuffer", { sharedArrayBuffer: false }, "the browser has no SharedArrayBuffer"],
    [
      "a browser that cannot create shared memory",
      { sharedMemoryFails: true },
      "the browser cannot create shared WebAssembly memory: Out of memory: Cannot allocate Wasm memory for new instance",
    ],
  ])("%s selects the single-threaded module with a warning", (_name, overrides, reason) => {
    expect(chooseRunThreads(8, scope(overrides), THREADED_WASM_URL)).toStrictEqual({
      jobs: 8,
      threadedWasmUrl: undefined,
      warning: `runs use one thread, because ${reason}`,
    });
  });
});

function scope({
  crossOriginIsolated = true,
  sharedArrayBuffer = true,
  sharedMemoryFails = false,
}: {
  crossOriginIsolated?: boolean;
  sharedArrayBuffer?: boolean;
  sharedMemoryFails?: boolean;
} = {}): ThreadScope {
  return {
    crossOriginIsolated,
    SharedArrayBuffer: sharedArrayBuffer ? SharedArrayBuffer : undefined,
    WebAssembly: {
      Memory: class {
        readonly buffer: ArrayBufferLike;

        constructor(descriptor: WebAssembly.MemoryDescriptor) {
          if (sharedMemoryFails || descriptor.maximum !== 65_536) {
            throw new RangeError("Out of memory: Cannot allocate Wasm memory for new instance");
          }

          this.buffer = descriptor.shared === true ? new SharedArrayBuffer(8) : new ArrayBuffer(8);
        }
      },
    },
  };
}
