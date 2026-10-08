import { errorMessage } from "@neherlab/app-contracts";
import createBackendModule from "@neherlab/app-wasm";

import { instantiateBackend } from "./instantiate";
import { type WorkerEvent, workerStartSchema } from "./protocol";
import { ExchangeServer, messagePortEndpoint, type WasmBackend } from "./serve";

const NO_PANIC_TEXT_MESSAGE =
  "the back end stopped without an error message; it most likely ran out of memory, because the browser limits the memory of a web page";

const scope: WorkerScope = self;

scope.addEventListener(
  "message",
  (event) => {
    void start(event.data);
  },
  { once: true },
);

// oxlint-disable-next-line anti-slop/no-unknown-parameters -- the first message of the page arrives untyped; it is parsed with the start schema
async function start(data: unknown): Promise<void> {
  const { module, port } = workerStartSchema.parse(data);
  let panicText: string | undefined;
  let backend: WasmBackend;

  try {
    const wasm = await instantiateBackend(module, createBackendModule);

    wasm.setPanicSink((text) => {
      panicText = text;
    });
    backend = new wasm.Backend();
  } catch (error) {
    post({ kind: "failed", message: errorMessage(error) });

    return;
  }

  const server = new ExchangeServer(messagePortEndpoint(port), backend, () => {
    post({ kind: "crashed", message: panicText ?? NO_PANIC_TEXT_MESSAGE });
  });

  server.listen();
  post({ kind: "ready" });
}

function post(event: WorkerEvent): void {
  // oxlint-disable-next-line unicorn/require-post-message-target-origin -- the scope of a dedicated worker takes no target origin
  scope.postMessage(event);
}

interface WorkerScope {
  addEventListener(type: "message", listener: (event: { data: unknown }) => void, options: { once: boolean }): void;
  postMessage(event: WorkerEvent): void;
}
