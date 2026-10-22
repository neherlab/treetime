import pWaitFor from "p-wait-for";
import { match } from "ts-pattern";

import type { EmscriptenOptions } from "../../backend/instantiate";
import type { JobMessage, ModuleVariant } from "../protocol";
import {
  type JobModule,
  type JobModuleFactory,
  type JobScope,
  serveJob,
  type WasmChunkReader,
  type WasmJob,
} from "../serve";

const EMPTY_MODULE = new WebAssembly.Module(new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]));

const START = { kind: "start", attempt: 3, variant: "single", module: EMPTY_MODULE, job: '{"id":"r1"}', jobs: 1 };

describe("job worker protocol", () => {
  test("a started job writes its inputs, streams its events, and sends its outputs before done", async () => {
    const { scope, runtime } = served({ outputs: { "ancestral.nwk": ["(A,", "B);"], "empty.txt": [] } });

    scope.receive(START);
    await scope.received("started");
    scope.receive({ kind: "input-begin", attempt: 3, path: "/in/tree.nwk" });
    scope.receive({ kind: "input-chunk", attempt: 3, path: "/in/tree.nwk", data: new TextEncoder().encode("(A,B);") });
    scope.receive({ kind: "run", attempt: 3 });
    await scope.received("done");

    expect([runtime.inputs, scope.posted.map(describeMessage)]).toStrictEqual([
      { "/in/tree.nwk": "(A,B);" },
      [
        "started",
        "input-ack",
        'event {"type":"log"}',
        "output-begin ancestral.nwk",
        "output-chunk ancestral.nwk (A,",
        "output-chunk ancestral.nwk B);",
        "output-begin empty.txt",
        "done {}",
      ],
    ]);
  });

  test("the worker holds at most 4 output chunks without an acknowledgement", async () => {
    const { scope } = served({ outputs: { "tree.nwk": ["1", "2", "3", "4", "5", "6"] } });

    scope.receive(START);
    await scope.received("started");
    scope.receive({ kind: "run", attempt: 3 });
    await pWaitFor(() => countKind(scope.posted, "output-chunk") === 4);
    scope.receive({ kind: "output-ack", attempt: 3 });
    await pWaitFor(() => countKind(scope.posted, "output-chunk") === 5);

    expect(scope.posted.map(({ kind }) => kind).includes("done")).toStrictEqual(false);
  });
});

describe("job worker failures", () => {
  test("an abort without stderr lines reports the abort reason as the crash", async () => {
    const { scope, runtime } = served();

    scope.receive(START);
    await scope.received("started");
    runtime.options()?.onAbort("stack overflow");

    expect(scope.posted.at(-1)).toStrictEqual({ kind: "crashed", attempt: 3, lines: ["stack overflow"] });
  });

  test("an abort after unrelated stderr lines keeps the stderr lines and the abort reason", async () => {
    const { scope, runtime } = served();

    scope.receive(START);
    await scope.received("started");
    runtime.options()?.printErr("warning: the tree has repeated node names");
    runtime.options()?.printErr("memory allocation of 4294967296 bytes failed");
    runtime.options()?.onAbort("Aborted()");

    expect(scope.posted.at(-1)).toStrictEqual({
      kind: "crashed",
      attempt: 3,
      lines: ["warning: the tree has repeated node names", "memory allocation of 4294967296 bytes failed", "Aborted()"],
    });
  });

  test("an error event without stderr lines reports its message and stays in the worker", async () => {
    const { scope } = served();

    scope.receive(START);
    await scope.received("started");
    const prevented = scope.error("RuntimeError: unreachable");

    expect([scope.posted.at(-1), prevented]).toStrictEqual([
      { kind: "crashed", attempt: 3, lines: ["RuntimeError: unreachable"] },
      true,
    ]);
  });

  test("a failure before the job started is a start failure", async () => {
    const { scope } = served({
      load: () => Promise.reject(new Error("this build has no threaded WebAssembly module")),
    });

    scope.receive({ ...START, variant: "threaded" });
    await scope.received("failed");

    expect(scope.posted).toStrictEqual([
      { kind: "failed", attempt: 3, lines: ["this build has no threaded WebAssembly module"] },
    ]);
  });

  test("a crash is reported once although the abort also throws out of the run", async () => {
    const { scope } = served({ abortInRun: "stack overflow" });

    scope.receive(START);
    await scope.received("started");
    scope.receive({ kind: "run", attempt: 3 });

    expect(scope.posted.filter(({ kind }) => kind === "crashed")).toStrictEqual([
      { kind: "crashed", attempt: 3, lines: ["stack overflow"] },
    ]);
  });
});

interface FakeScope extends JobScope {
  posted: JobMessage[];
  receive(data: unknown): void;
  error(message: string): boolean;
  received(kind: JobMessage["kind"]): Promise<void>;
}

interface ServedJob {
  scope: FakeScope;
  runtime: FakeRuntime;
}

interface FakeRuntime {
  inputs: Record<string, string>;
  options(): EmscriptenOptions | undefined;
}

function served({
  outputs = {},
  load,
  abortInRun,
}: {
  outputs?: Record<string, string[]>;
  load?: (variant: ModuleVariant) => Promise<JobModuleFactory>;
  abortInRun?: string;
} = {}): ServedJob {
  const scope = fakeScope();
  const inputs: Record<string, string> = {};
  let options: EmscriptenOptions | undefined;

  const job: WasmJob = {
    beginInput: (path) => {
      inputs[path] = "";
    },
    appendInput: (path, chunk) => {
      inputs[path] = `${inputs[path] ?? ""}${new TextDecoder().decode(chunk)}`;
    },
    run: (onEvent) => {
      if (abortInRun !== undefined) {
        options?.onAbort(abortInRun);

        throw new WebAssembly.RuntimeError(`Aborted(${abortInRun})`);
      }

      onEvent('{"type":"log"}');

      return Promise.resolve("{}");
    },
    outputFiles: () => Object.keys(outputs),
    readOutput: (path) => {
      const chunks = (outputs[path] ?? []).map((text) => new TextEncoder().encode(text));

      return { chunk: () => chunks.shift(), free: () => undefined };
    },
  };

  const module: JobModule = {
    Job: class implements WasmJob {
      beginInput(path: string): void {
        job.beginInput(path);
      }

      appendInput(path: string, chunk: Uint8Array): void {
        job.appendInput(path, chunk);
      }

      run(onEvent: (event: string) => void): Promise<string> {
        return job.run(onEvent);
      }

      outputFiles(): string[] {
        return job.outputFiles();
      }

      readOutput(path: string): WasmChunkReader {
        return job.readOutput(path);
      }
    },
  };

  const create: JobModuleFactory = (received) => {
    options = received;

    return Promise.resolve(module);
  };

  serveJob(scope, load ?? (() => Promise.resolve(create)));

  return { scope, runtime: { inputs, options: () => options } };
}

function fakeScope(): FakeScope {
  const messageListeners: Array<(event: { data: unknown }) => void> = [];
  const errorListeners: Array<(event: { message?: string | undefined; preventDefault(): void }) => void> = [];

  const scope: FakeScope = {
    posted: [],
    addEventListener(type: "message" | "error", listener: never) {
      if (type === "message") {
        messageListeners.push(listener);
      } else {
        errorListeners.push(listener);
      }
    },
    postMessage(message) {
      scope.posted.push(message);
    },
    receive(data) {
      messageListeners.forEach((listener) => {
        listener({ data });
      });
    },
    error(message) {
      let prevented = false;

      errorListeners.forEach((listener) => {
        listener({
          message,
          preventDefault: () => {
            prevented = true;
          },
        });
      });

      return prevented;
    },
    received: (kind) => pWaitFor(() => scope.posted.some((message) => message.kind === kind)),
  };

  return scope;
}

function countKind(messages: JobMessage[], kind: JobMessage["kind"]): number {
  return messages.filter((message) => message.kind === kind).length;
}

function describeMessage(message: JobMessage): string {
  return match(message)
    .with({ kind: "event" }, ({ event }) => `event ${event}`)
    .with({ kind: "output-begin" }, ({ path }) => `output-begin ${path}`)
    .with({ kind: "output-chunk" }, ({ path, data }) => `output-chunk ${path} ${new TextDecoder().decode(data)}`)
    .with({ kind: "done" }, ({ outcome }) => `done ${outcome}`)
    .otherwise(({ kind }) => kind);
}
