import { errorMessage } from "@neherlab/app-contracts";

import { type EmscriptenOptions, ErrorLines, instantiateBackend } from "./instantiate";
import { type JobBackend, JobHost, type JobModules, type JobWorker } from "./jobs";
import { type WorkerEvent, type WorkerStart, workerStartSchema } from "./protocol";
import { ExchangeServer, messagePortEndpoint, type WasmBackend } from "./serve";

export const NO_ERROR_LINES_MESSAGE = "the back end stopped without an error message";

export function serveBackend(
  scope: BackendScope,
  create: (options: EmscriptenOptions) => Promise<BackendModule>,
  startJobWorker: () => JobWorker,
  compile: (url: string) => Promise<WebAssembly.Module>,
): void {
  const backend = new BackendWorker(scope, create, startJobWorker, compile);

  scope.addEventListener(
    "message",
    (event) => {
      void backend.start(event.data);
    },
    { once: true },
  );
  scope.addEventListener("error", (event) => {
    event.preventDefault();
    backend.fail(event.message ?? "the back end worker stopped with an error");
  });
}

export interface BackendScope {
  addEventListener(type: "message", listener: (event: { data: unknown }) => void, options: { once: boolean }): void;
  addEventListener(
    type: "error",
    listener: (event: { message?: string | undefined; preventDefault(): void }) => void,
  ): void;
  postMessage(event: WorkerEvent): void;
}

export interface BackendModule {
  Backend: new () => WasmBackend & JobBackend;
}

class BackendWorker {
  readonly #scope: BackendScope;
  readonly #create: (options: EmscriptenOptions) => Promise<BackendModule>;
  readonly #startJobWorker: () => JobWorker;
  readonly #compile: (url: string) => Promise<WebAssembly.Module>;
  readonly #errors = new ErrorLines(() => {
    this.#end();
  });
  #ready = false;
  #ended = false;

  constructor(
    scope: BackendScope,
    create: (options: EmscriptenOptions) => Promise<BackendModule>,
    startJobWorker: () => JobWorker,
    compile: (url: string) => Promise<WebAssembly.Module>,
  ) {
    this.#scope = scope;
    this.#create = create;
    this.#startJobWorker = startJobWorker;
    this.#compile = compile;
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- the first message of the page arrives untyped; it is parsed with the start schema
  async start(data: unknown): Promise<void> {
    let backend: WasmBackend & JobBackend;
    let port: MessagePort;

    try {
      const start = workerStartSchema.parse(data);
      const wasm = await instantiateBackend(start.module, this.#create, this.#errors);

      backend = new wasm.Backend();
      port = start.port;

      const onTrap = (trap: WebAssembly.RuntimeError): void => {
        this.fail(trap.message);
      };

      const jobs = new JobHost({
        backend,
        startWorker: this.#startJobWorker,
        modules: this.#jobModules(start),
        onTrap,
      });

      const server = new ExchangeServer(messagePortEndpoint(port), backend, onTrap, () => {
        jobs.drain();
      });

      server.listen();
    } catch (error) {
      this.fail(errorMessage(error));

      return;
    }

    if (!this.#ended) {
      this.#ready = true;
      this.#post({ kind: "ready" });
    }
  }

  fail(message: string): void {
    this.#errors.error(message);
    this.#end();
  }

  #jobModules({ module, jobs, threadedWasmUrl }: WorkerStart): JobModules {
    if (threadedWasmUrl === undefined) {
      return { single: module };
    }

    return { single: module, threaded: { module: () => this.#compile(threadedWasmUrl), jobs } };
  }

  #end(): void {
    if (this.#ended) {
      return;
    }

    this.#ended = true;

    const lines = this.#errors.lines();
    const message = lines.length === 0 ? NO_ERROR_LINES_MESSAGE : lines.join("\n");

    this.#post({ kind: this.#ready ? "crashed" : "failed", message });
  }

  #post(event: WorkerEvent): void {
    // oxlint-disable-next-line unicorn/require-post-message-target-origin -- the scope of a dedicated worker takes no target origin
    this.#scope.postMessage(event);
  }
}
