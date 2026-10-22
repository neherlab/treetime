import { errorMessage } from "@neherlab/app-contracts";
import type { BackendStopped } from "@neherlab/app-ui/host";
import type { FetchConnection, FetchPort } from "@neherlab/app-ui/port-fetch";
import { match } from "ts-pattern";

import { type WorkerEvent, type WorkerStart, workerEventSchema } from "./protocol";

export class WorkerConnection implements FetchConnection {
  readonly #host: WorkerHost;
  readonly #portListeners = new Set<(port: FetchPort) => void>();
  readonly #stopListeners = new Set<(stop: BackendStopped) => void>();
  #module: Promise<WebAssembly.Module> | undefined;
  #worker: RunningWorker | undefined;
  #stopped = false;

  constructor(host: WorkerHost) {
    this.#host = host;
  }

  onPort(listener: (port: FetchPort) => void): void {
    this.#portListeners.add(listener);

    const port = this.#worker?.port;

    if (port === undefined) {
      this.#start();
    } else {
      listener(port);
    }
  }

  onStopped(listener: (stop: BackendStopped) => void): void {
    this.#stopListeners.add(listener);
  }

  #start(): void {
    if (this.#worker !== undefined || this.#stopped) {
      return;
    }

    const channel = this.#host.channel();
    const running: RunningWorker = { worker: this.#host.start(), port: undefined };

    this.#worker = running;
    running.worker.addEventListener("message", (event) => {
      this.#receive(running, channel.port1, event.data);
    });
    running.worker.addEventListener("error", (event) => {
      this.#crash(running, event.message ?? "the back end worker stopped with an error");
    });
    void this.#send(running, channel.port2);
  }

  async #send(running: RunningWorker, port: MessagePort): Promise<void> {
    let module: WebAssembly.Module;

    try {
      module = await this.#compiled();
    } catch (error) {
      this.#fail(running, `loading the WebAssembly module failed: ${errorMessage(error)}`);

      return;
    }

    if (this.#worker === running) {
      const { jobs, threadedWasmUrl } = this.#host.runThreads;

      running.worker.postMessage({ kind: "start", module, port, jobs, threadedWasmUrl }, [port]);
    }
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- worker messages arrive untyped; they are parsed with the event schema
  #receive(running: RunningWorker, port: FetchPort, data: unknown): void {
    const event = workerEventSchema.safeParse(data);

    if (!event.success) {
      this.#fail(running, "the back end worker sent a malformed message");

      return;
    }

    match<WorkerEvent>(event.data)
      .with({ kind: "ready" }, () => {
        this.#ready(running, port);
      })
      .with({ kind: "failed" }, ({ message }) => {
        this.#fail(running, message);
      })
      .with({ kind: "crashed" }, ({ message }) => {
        this.#crash(running, message);
      })
      .exhaustive();
  }

  #ready(running: RunningWorker, port: FetchPort): void {
    if (this.#worker !== running) {
      return;
    }

    running.port = port;
    this.#portListeners.forEach((listener) => {
      listener(port);
    });
  }

  #crash(running: RunningWorker, message: string): void {
    if (this.#worker !== running) {
      return;
    }

    if (running.port === undefined) {
      this.#fail(running, message);

      return;
    }

    this.#end(running, { reason: `the back end stopped and restarts: ${message}`, restarts: true });
    this.#start();
  }

  #fail(running: RunningWorker, message: string): void {
    if (this.#worker !== running) {
      return;
    }

    this.#stopped = true;
    this.#end(running, { reason: `the back end could not start: ${message}`, restarts: false });
  }

  #end(running: RunningWorker, stop: BackendStopped): void {
    this.#worker = undefined;
    running.worker.terminate();
    this.#stopListeners.forEach((listener) => {
      listener(stop);
    });
  }

  #compiled(): Promise<WebAssembly.Module> {
    if (this.#module === undefined) {
      const module = this.#host.compile();

      this.#module = module;
      module.catch(() => {
        if (this.#module === module) {
          this.#module = undefined;
        }
      });
    }

    return this.#module;
  }
}

export interface WorkerHost {
  start(): BackendWorker;
  compile(): Promise<WebAssembly.Module>;
  channel(): WorkerChannel;
  readonly runThreads: WorkerRunThreads;
}

export interface WorkerRunThreads {
  jobs: number;
  threadedWasmUrl: string | undefined;
}

export interface BackendWorker {
  postMessage(message: WorkerStart, transfer: Transferable[]): void;
  addEventListener(type: "message", listener: (event: { data: unknown }) => void): void;
  addEventListener(type: "error", listener: (event: { message?: string | undefined }) => void): void;
  terminate(): void;
}

export interface WorkerChannel {
  port1: FetchPort;
  port2: MessagePort;
}

interface RunningWorker {
  worker: BackendWorker;
  port: FetchPort | undefined;
}
