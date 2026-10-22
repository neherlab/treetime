import pWaitFor from "p-wait-for";

import type { HostMessage } from "../../job/protocol";
import { type BackendModule, type BackendScope, NO_ERROR_LINES_MESSAGE, serveBackend } from "../backend-worker";
import type { EmscriptenOptions } from "../instantiate";
import type { JobWorker, WasmJobOrder } from "../jobs";
import type { WorkerEvent } from "../protocol";
import type { WasmExchange } from "../serve";

const EMPTY_MODULE = new WebAssembly.Module(new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]));

const THREADED_MODULE = new WebAssembly.Module(new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]));

const THREADED_WASM_URL = "https://treetime.localhost/assets/app_wasm-mt.wasm";

describe("backend worker run modules", () => {
  test("runs start on the threaded module of the start message with its thread count", async () => {
    const { posted, compiled, close } = servedRun({ jobs: 8, threadedWasmUrl: THREADED_WASM_URL });

    await pWaitFor(() => posted.length === 1);
    close();

    expect([compiled, posted]).toStrictEqual([
      [THREADED_WASM_URL],
      [{ kind: "start", attempt: 1, variant: "threaded", module: THREADED_MODULE, job: "{}", jobs: 8 }],
    ]);
  });

  test("runs start on the single-threaded module with one thread when the start message has no threaded module", async () => {
    const { posted, compiled, close } = servedRun({ jobs: 8 });

    await pWaitFor(() => posted.length === 1);
    close();

    expect([compiled, posted]).toStrictEqual([
      [],
      [{ kind: "start", attempt: 1, variant: "single", module: EMPTY_MODULE, job: "{}", jobs: 1 }],
    ]);
  });
});

describe("backend worker failures", () => {
  test("an abort after the start reports the stderr lines and the abort reason as the crash message", async () => {
    const { scope, options, close } = served();

    await pWaitFor(() => scope.posted.length === 1);
    options()?.printErr("thread 'main' panicked at src/run.rs:1:1: index out of bounds");
    options()?.onAbort("Aborted()");
    close();

    expect(scope.posted).toStrictEqual([
      { kind: "ready" },
      { kind: "crashed", message: "thread 'main' panicked at src/run.rs:1:1: index out of bounds\nAborted()" },
    ]);
  });

  test("an abort without any error text reports that the back end stopped without a message", async () => {
    const { scope, options, close } = served();

    await pWaitFor(() => scope.posted.length === 1);
    options()?.onAbort(undefined);
    close();

    expect(scope.posted.at(-1)).toStrictEqual({ kind: "crashed", message: NO_ERROR_LINES_MESSAGE });
  });

  test("a module that cannot start is a failed start with its error", async () => {
    const { scope, close } = served(() => Promise.reject(new Error("When creating the runs directory: unsupported")));

    await pWaitFor(() => scope.posted.length === 1);
    close();

    expect(scope.posted).toStrictEqual([{ kind: "failed", message: "When creating the runs directory: unsupported" }]);
  });
});

interface FakeScope extends BackendScope {
  posted: WorkerEvent[];
}

interface ServedBackend {
  scope: FakeScope;
  options: () => EmscriptenOptions | undefined;
  close: () => void;
}

function served(create?: (options: EmscriptenOptions) => Promise<BackendModule>): ServedBackend {
  const channel = new MessageChannel();
  let options: EmscriptenOptions | undefined;
  let start: (event: { data: unknown }) => void = () => undefined;

  const scope: FakeScope = {
    posted: [],
    addEventListener(type: "message" | "error", listener: never) {
      if (type === "message") {
        start = listener;
      }
    },
    postMessage(event) {
      scope.posted.push(event);
    },
  };

  const createModule = (received: EmscriptenOptions): Promise<BackendModule> => {
    options = received;

    return Promise.resolve({ Backend: FakeBackend });
  };

  serveBackend(scope, create ?? createModule, unusedJobWorker, unusedCompile);
  start({ data: { kind: "start", module: EMPTY_MODULE, port: channel.port2, jobs: 1 } });

  return {
    scope,
    options: () => options,
    close: () => {
      channel.port1.close();
      channel.port2.close();
    },
  };
}

function unusedJobWorker(): JobWorker {
  throw new Error("the test starts no run");
}

function unusedCompile(): Promise<WebAssembly.Module> {
  throw new Error("the test compiles no threaded module");
}

interface ServedRun {
  posted: HostMessage[];
  compiled: string[];
  close: () => void;
}

function servedRun(threads: { jobs: number; threadedWasmUrl?: string }): ServedRun {
  const channel = new MessageChannel();
  const posted: HostMessage[] = [];
  const compiled: string[] = [];
  let start: (event: { data: unknown }) => void = () => undefined;

  const scope: BackendScope = {
    addEventListener(type: "message" | "error", listener: never) {
      if (type === "message") {
        start = listener;
      }
    },
    postMessage: () => undefined,
  };

  const jobWorker: JobWorker = {
    postMessage: (message) => {
      posted.push(message);
    },
    addEventListener: () => undefined,
    terminate: () => undefined,
  };

  serveBackend(
    scope,
    () => Promise.resolve({ Backend: RunningBackend }),
    () => jobWorker,
    (url) => {
      compiled.push(url);

      return Promise.resolve(THREADED_MODULE);
    },
  );
  start({ data: { kind: "start", module: EMPTY_MODULE, port: channel.port2, ...threads } });
  channel.port1.postMessage({
    kind: "request",
    request: { seq: 0, method: "POST", url: "/api/runs/r1/start", headers: [] },
  });

  return {
    posted,
    compiled,
    close: () => {
      channel.port1.close();
      channel.port2.close();
    },
  };
}

class FakeBackend {
  fetch(): WasmExchange {
    throw new Error("the test sends no request");
  }

  takeJobOrders(): WasmJobOrder[] {
    return [];
  }

  readRunInput(): never {
    throw new Error("the test reads no input");
  }

  appendRunEvent(): void {
    throw new Error("the test records no event");
  }

  beginRunOutput(): void {
    throw new Error("the test writes no output");
  }

  appendRunOutput(): void {
    throw new Error("the test writes no output");
  }

  completeRun(): boolean {
    throw new Error("the test completes no run");
  }

  failRun(): boolean {
    throw new Error("the test fails no run");
  }
}

class RunningBackend extends FakeBackend {
  #orders = [{ kind: "start", runId: "r1", job: "{}", inputPaths: [], free: () => undefined }];

  override fetch(): WasmExchange {
    return {
      head: () => Promise.resolve({ status: 200, headerNames: [], headerValues: [] }),
      chunk: () => Promise.resolve(undefined),
      abort: () => undefined,
      free: () => undefined,
    };
  }

  override takeJobOrders(): WasmJobOrder[] {
    return this.#orders.splice(0);
  }
}
