import type { BackendStopped } from "@neherlab/app-ui/host";
import type { FetchPort } from "@neherlab/app-ui/port-fetch";
import { describe, expect, test } from "vitest";

import { type BackendWorker, WorkerConnection, type WorkerHost } from "../connection";
import type { WorkerStart } from "../protocol";

const EMPTY_MODULE = new WebAssembly.Module(new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]));

describe("connection start", () => {
  test("the first port listener starts a worker with the compiled module and the worker end of a channel", async () => {
    const host = fakeHost();
    const connection = new WorkerConnection(host);

    connection.onPort(() => undefined);
    const [worker] = host.workers;
    await worker?.started();

    expect([host.workers.length, worker?.posted, worker?.transferred]).toStrictEqual([
      1,
      [{ kind: "start", module: EMPTY_MODULE, port: host.channels[0]?.port2 }],
      [[host.channels[0]?.port2]],
    ]);
  });

  test("a ready worker gives every port listener the page end of the channel", async () => {
    const host = fakeHost();
    const connection = new WorkerConnection(host);
    const ports: FetchPort[] = [];

    connection.onPort((port) => {
      ports.push(port);
    });
    await host.workers[0]?.started();
    host.workers[0]?.emit({ kind: "ready" });
    connection.onPort((port) => {
      ports.push(port);
    });

    expect(ports).toStrictEqual([host.channels[0]?.port1, host.channels[0]?.port1]);
  });
});

describe("connection failures", () => {
  test("a worker that cannot start stops the back end for good", async () => {
    const host = fakeHost();
    const connection = new WorkerConnection(host);
    const stops: BackendStopped[] = [];

    connection.onStopped((stop) => {
      stops.push(stop);
    });
    connection.onPort(() => undefined);
    await host.workers[0]?.started();
    host.workers[0]?.emit({
      kind: "failed",
      message: "When creating the runs directory '/treetime/runs': unsupported",
    });
    connection.onPort(() => undefined);

    expect([stops, host.workers.length, host.workers[0]?.terminated]).toStrictEqual([
      [
        {
          reason: "the back end could not start: When creating the runs directory '/treetime/runs': unsupported",
          restarts: false,
        },
      ],
      1,
      true,
    ]);
  });

  test("a module that does not load stops the back end for good", async () => {
    const host = fakeHost({ compile: () => Promise.reject(new Error("HTTP status 404")) });
    const connection = new WorkerConnection(host);
    const stops: BackendStopped[] = [];

    connection.onStopped((stop) => {
      stops.push(stop);
    });
    connection.onPort(() => undefined);
    await expect.poll(() => stops).toHaveLength(1);

    expect([stops, host.workers[0]?.posted]).toStrictEqual([
      [
        {
          reason: "the back end could not start: loading the WebAssembly module failed: HTTP status 404",
          restarts: false,
        },
      ],
      [],
    ]);
  });

  test("a crash after the start restarts the worker with the same module", async () => {
    const host = fakeHost();
    const connection = new WorkerConnection(host);
    const stops: BackendStopped[] = [];
    const ports: FetchPort[] = [];

    connection.onStopped((stop) => {
      stops.push(stop);
    });
    connection.onPort((port) => {
      ports.push(port);
    });
    await host.workers[0]?.started();
    host.workers[0]?.emit({ kind: "ready" });
    host.workers[0]?.emit({ kind: "crashed", message: "panicked at src/run.rs:1:1: index out of bounds" });
    await host.workers[1]?.started();
    host.workers[1]?.emit({ kind: "ready" });

    expect([stops, host.compiled, host.workers[0]?.terminated, ports]).toStrictEqual([
      [
        {
          reason: "the back end stopped and restarts: panicked at src/run.rs:1:1: index out of bounds",
          restarts: true,
        },
      ],
      1,
      true,
      [host.channels[0]?.port1, host.channels[1]?.port1],
    ]);
  });

  test("an error of a running worker restarts it", async () => {
    const host = fakeHost();
    const connection = new WorkerConnection(host);
    const stops: BackendStopped[] = [];

    connection.onStopped((stop) => {
      stops.push(stop);
    });
    connection.onPort(() => undefined);
    await host.workers[0]?.started();
    host.workers[0]?.emit({ kind: "ready" });
    host.workers[0]?.fail("Uncaught ZodError");

    expect([stops, host.workers.length]).toStrictEqual([
      [{ reason: "the back end stopped and restarts: Uncaught ZodError", restarts: true }],
      2,
    ]);
  });

  test("a crash before the worker is ready is a failed start", async () => {
    const host = fakeHost();
    const connection = new WorkerConnection(host);
    const stops: BackendStopped[] = [];

    connection.onStopped((stop) => {
      stops.push(stop);
    });
    connection.onPort(() => undefined);
    await host.workers[0]?.started();
    host.workers[0]?.fail("SyntaxError: Unexpected token");

    expect([stops, host.workers.length]).toStrictEqual([
      [{ reason: "the back end could not start: SyntaxError: Unexpected token", restarts: false }],
      1,
    ]);
  });

  test("a malformed message of the worker is a failed start", async () => {
    const host = fakeHost();
    const connection = new WorkerConnection(host);
    const stops: BackendStopped[] = [];

    connection.onStopped((stop) => {
      stops.push(stop);
    });
    connection.onPort(() => undefined);
    await host.workers[0]?.started();
    host.workers[0]?.emit({ kind: "ready", extra: true });

    expect(stops).toStrictEqual([
      { reason: "the back end could not start: the back end worker sent a malformed message", restarts: false },
    ]);
  });
});

interface FakeHost extends WorkerHost {
  workers: FakeWorker[];
  channels: MessageChannel[];
  readonly compiled: number;
}

interface FakeWorker extends BackendWorker {
  posted: WorkerStart[];
  transferred: Transferable[][];
  terminated: boolean;
  started(): Promise<void>;
  emit(data: unknown): void;
  fail(message: string): void;
}

function fakeHost({ compile }: { compile?: () => Promise<WebAssembly.Module> } = {}): FakeHost {
  let compiled = 0;

  const host: FakeHost = {
    workers: [],
    channels: [],
    get compiled() {
      return compiled;
    },
    start: () => {
      const worker = fakeWorker();

      host.workers.push(worker);

      return worker;
    },
    compile: () => {
      compiled += 1;

      return compile === undefined ? Promise.resolve(EMPTY_MODULE) : compile();
    },
    channel: () => {
      const channel = new MessageChannel();

      host.channels.push(channel);

      return channel;
    },
  };

  return host;
}

function fakeWorker(): FakeWorker {
  const messageListeners: Array<(event: { data: unknown }) => void> = [];
  const errorListeners: Array<(event: { message?: string | undefined }) => void> = [];
  const { promise: startPosted, resolve: postStart } = Promise.withResolvers<undefined>();

  const worker: FakeWorker = {
    posted: [],
    transferred: [],
    terminated: false,
    postMessage(message, transfer) {
      worker.posted.push(message);
      worker.transferred.push(transfer);
      postStart(undefined);
    },
    addEventListener(type: "message" | "error", listener: never) {
      if (type === "message") {
        messageListeners.push(listener);
      } else {
        errorListeners.push(listener);
      }
    },
    terminate() {
      worker.terminated = true;
    },
    started: () => startPosted,
    emit(data) {
      messageListeners.forEach((listener) => {
        listener({ data });
      });
    },
    fail(message) {
      errorListeners.forEach((listener) => {
        listener({ message });
      });
    },
  };

  return worker;
}
