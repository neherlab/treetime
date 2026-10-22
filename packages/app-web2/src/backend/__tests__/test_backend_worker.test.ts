import pWaitFor from "p-wait-for";

import { type BackendModule, type BackendScope, NO_ERROR_LINES_MESSAGE, serveBackend } from "../backend-worker";
import type { EmscriptenOptions } from "../instantiate";
import type { JobWorker } from "../jobs";
import type { WorkerEvent } from "../protocol";

const EMPTY_MODULE = new WebAssembly.Module(new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]));

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

  serveBackend(scope, create ?? createModule, unusedJobWorker);
  start({ data: { kind: "start", module: EMPTY_MODULE, port: channel.port2 } });

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

class FakeBackend {
  fetch(): never {
    throw new Error("the test sends no request");
  }

  takeJobOrders(): [] {
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
