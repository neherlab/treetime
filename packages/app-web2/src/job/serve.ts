import { errorMessage } from "@neherlab/app-contracts";
import { match } from "ts-pattern";

import { type EmscriptenOptions, ErrorLines, instantiateBackend } from "../backend/instantiate";
import {
  ChunkWindow,
  type HostMessage,
  hostMessageSchema,
  type JobMessage,
  messageAttemptSchema,
  type ModuleVariant,
  transferList,
} from "./protocol";

export function serveJob(scope: JobScope, load: (variant: ModuleVariant) => Promise<JobModuleFactory>): void {
  const server = new JobServer(scope, load);

  scope.addEventListener("message", (event) => {
    server.receive(event.data);
  });
  scope.addEventListener("error", (event) => {
    event.preventDefault();
    server.fail(event.message ?? "the run worker stopped with an error");
  });
}

export interface JobScope {
  addEventListener(type: "message", listener: (event: { data: unknown }) => void): void;
  addEventListener(
    type: "error",
    listener: (event: { message?: string | undefined; preventDefault(): void }) => void,
  ): void;
  postMessage(message: JobMessage, options: { transfer: Transferable[] }): void;
}

export type JobModuleFactory = (options: EmscriptenOptions) => Promise<JobModule>;

export interface JobModule {
  Job: new (jobJson: string, jobs: number) => WasmJob;
}

export interface WasmJob {
  beginInput(path: string): void;
  appendInput(path: string, chunk: Uint8Array): void;
  run(onEvent: (event: string) => void): Promise<string>;
  outputFiles(): string[];
  readOutput(path: string): WasmChunkReader;
}

export interface WasmChunkReader {
  chunk(): Uint8Array | undefined;
  free(): void;
}

class JobServer {
  readonly #scope: JobScope;
  readonly #load: (variant: ModuleVariant) => Promise<JobModuleFactory>;
  readonly #errors = new ErrorLines(() => {
    this.#end();
  });
  readonly #window = new ChunkWindow();
  #attempt: number | undefined;
  #job: WasmJob | undefined;
  #started = false;
  #ended = false;

  constructor(scope: JobScope, load: (variant: ModuleVariant) => Promise<JobModuleFactory>) {
    this.#scope = scope;
    this.#load = load;
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- messages of the job host arrive untyped; they are parsed with the host message schema
  receive(data: unknown): void {
    this.#attempt ??= messageAttemptSchema.safeParse(data).data?.attempt;

    try {
      this.#handle(hostMessageSchema.parse(data));
    } catch (error) {
      this.fail(errorMessage(error));
    }
  }

  fail(message: string): void {
    this.#errors.error(message);
    this.#end();
  }

  #handle(message: HostMessage): void {
    match(message)
      .with({ kind: "start" }, (start) => {
        void this.#start(start);
      })
      .with({ kind: "input-begin" }, ({ path }) => {
        this.#running().beginInput(path);
      })
      .with({ kind: "input-chunk" }, ({ path, data }) => {
        this.#running().appendInput(path, data);
        this.#post({ kind: "input-ack", attempt: message.attempt });
      })
      .with({ kind: "run" }, () => {
        void this.#run(this.#running());
      })
      .with({ kind: "output-ack" }, () => {
        this.#window.release();
      })
      .exhaustive();
  }

  async #start({ variant, module, job, jobs }: Extract<HostMessage, { kind: "start" }>): Promise<void> {
    try {
      const create = await this.#load(variant);
      const wasm = await instantiateBackend(module, create, this.#errors);

      this.#job = new wasm.Job(job, jobs);
    } catch (error) {
      this.fail(errorMessage(error));

      return;
    }

    if (!this.#ended) {
      this.#started = true;
      this.#post({ kind: "started", attempt: this.#current() });
    }
  }

  async #run(job: WasmJob): Promise<void> {
    try {
      const attempt = this.#current();

      const outcome = await job.run((event) => {
        this.#post({ kind: "event", attempt, event });
      });

      for (const path of job.outputFiles()) {
        if (!(await this.#sendOutput(job, attempt, path))) {
          return;
        }
      }

      this.#ended = true;
      this.#post({ kind: "done", attempt, outcome });
    } catch (error) {
      this.fail(errorMessage(error));
    }
  }

  async #sendOutput(job: WasmJob, attempt: number, path: string): Promise<boolean> {
    const reader = job.readOutput(path);

    try {
      this.#post({ kind: "output-begin", attempt, path });

      for (let chunk = reader.chunk(); chunk !== undefined; chunk = reader.chunk()) {
        if (!(await this.#window.acquire())) {
          return false;
        }

        this.#post({ kind: "output-chunk", attempt, path, data: chunk }, transferList(chunk));
      }

      return true;
    } finally {
      reader.free();
    }
  }

  #running(): WasmJob {
    if (this.#job === undefined) {
      throw new Error("the run worker received a run message before its module started");
    }

    return this.#job;
  }

  #current(): number {
    if (this.#attempt === undefined) {
      throw new Error("the run worker received no start message");
    }

    return this.#attempt;
  }

  #end(): void {
    if (this.#ended) {
      return;
    }

    this.#ended = true;
    this.#window.close();

    if (this.#attempt !== undefined) {
      const lines = this.#errors.lines();

      this.#post({ kind: this.#started ? "crashed" : "failed", attempt: this.#attempt, lines });
    }
  }

  #post(message: JobMessage, transfer: Transferable[] = []): void {
    this.#scope.postMessage(message, { transfer });
  }
}
