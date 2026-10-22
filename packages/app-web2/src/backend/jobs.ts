import { errorMessage } from "@neherlab/app-contracts";
import { match } from "ts-pattern";

import {
  ChunkWindow,
  type HostMessage,
  type JobMessage,
  jobMessageSchema,
  type ModuleVariant,
  transferList,
} from "../job/protocol";
import type { WasmChunkReader } from "../job/serve";

export const JOB_START_TIMEOUT_MS = 60_000;

const SINGLE_THREADED_JOBS = 1;

const NO_WORKER_ERROR_MESSAGE = "the run worker stopped with an error";

const MALFORMED_MESSAGE = "the run worker sent a malformed message";

export class JobHost {
  readonly #backend: JobBackend;
  readonly #startWorker: () => JobWorker;
  readonly #modules: JobModules;
  readonly #onTrap: (trap: WebAssembly.RuntimeError) => void;
  readonly #startTimeoutMs: number;
  readonly #jobs = new Map<string, RunJob>();
  #attempts = 0;
  #threadedFailed = false;
  #threadedModule: Promise<WebAssembly.Module> | undefined;
  #trapped = false;

  constructor(options: {
    backend: JobBackend;
    startWorker: () => JobWorker;
    modules: JobModules;
    onTrap: (trap: WebAssembly.RuntimeError) => void;
    startTimeoutMs?: number | undefined;
  }) {
    this.#backend = options.backend;
    this.#startWorker = options.startWorker;
    this.#modules = options.modules;
    this.#onTrap = options.onTrap;
    this.#startTimeoutMs = options.startTimeoutMs ?? JOB_START_TIMEOUT_MS;
  }

  drain(): void {
    if (this.#trapped) {
      return;
    }

    let orders: JobOrder[];

    try {
      orders = this.#backend.takeJobOrders().map(readOrder);
    } catch (error) {
      this.#backendFailed(error);

      return;
    }

    orders.forEach((order) => {
      match(order)
        .with({ kind: "start" }, ({ runId, job, inputPaths }) => {
          this.#begin(runId, job, inputPaths);
        })
        .with({ kind: "cancel" }, ({ runId }) => {
          this.#cancel(runId);
        })
        .exhaustive();
    });
  }

  #begin(runId: string, job: string, inputPaths: string[]): void {
    const run: RunJob = {
      runId,
      job,
      inputPaths,
      state: "starting",
      attempt: 0,
      variant: "single",
      worker: undefined,
      window: new ChunkWindow(),
      timer: undefined,
    };

    this.#jobs.set(runId, run);
    this.#launch(run);
  }

  #launch(run: RunJob): void {
    const attempt = this.#attempts + 1;
    const variant: ModuleVariant = this.#modules.threaded === undefined || this.#threadedFailed ? "single" : "threaded";

    this.#attempts = attempt;
    run.attempt = attempt;
    run.variant = variant;
    run.state = "starting";
    run.window = new ChunkWindow();

    let worker: JobWorker;

    try {
      worker = this.#startWorker();
    } catch (error) {
      this.#startFailed(run, attempt, [errorMessage(error)]);

      return;
    }

    run.worker = worker;
    worker.addEventListener("message", (event) => {
      this.#receive(run, attempt, event.data);
    });
    worker.addEventListener("error", (event) => {
      this.#lost(run, attempt, [event.message ?? NO_WORKER_ERROR_MESSAGE]);
    });
    run.timer = setTimeout(() => {
      this.#startFailed(run, attempt, [`the run worker did not start within ${String(this.#startTimeoutMs / 1000)} s`]);
    }, this.#startTimeoutMs);
    void this.#sendStart(run, attempt, variant);
  }

  async #sendStart(run: RunJob, attempt: number, variant: ModuleVariant): Promise<void> {
    let module: WebAssembly.Module;
    let jobs: number;

    try {
      [module, jobs] = await this.#module(variant);
    } catch (error) {
      this.#startFailed(run, attempt, [`loading the ${variant} WebAssembly module failed: ${errorMessage(error)}`]);

      return;
    }

    if (this.#current(run, attempt)) {
      this.#post(run, { kind: "start", attempt, variant, module, job: run.job, jobs });
    }
  }

  async #module(variant: ModuleVariant): Promise<[WebAssembly.Module, number]> {
    const { single, threaded } = this.#modules;

    if (variant === "single" || threaded === undefined) {
      return [single, SINGLE_THREADED_JOBS];
    }

    this.#threadedModule ??= threaded.module();

    return [await this.#threadedModule, threaded.jobs];
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- worker messages arrive untyped; they are parsed with the job message schema
  #receive(run: RunJob, attempt: number, data: unknown): void {
    if (!this.#current(run, attempt)) {
      return;
    }

    const parsed = jobMessageSchema.safeParse(data);

    if (!parsed.success) {
      this.#lost(run, attempt, [MALFORMED_MESSAGE]);

      return;
    }

    if (parsed.data.attempt === attempt) {
      this.#handle(run, attempt, parsed.data);
    }
  }

  #handle(run: RunJob, attempt: number, message: JobMessage): void {
    match(message)
      .with({ kind: "started" }, () => {
        if (run.state === "starting") {
          clearTimeout(run.timer);
          run.state = "inputs";
          void this.#sendInputs(run, attempt);
        }
      })
      .with({ kind: "failed" }, { kind: "crashed" }, ({ lines }) => {
        this.#lost(run, attempt, lines);
      })
      .with({ kind: "input-ack" }, () => {
        run.window.release();
      })
      .with({ kind: "event" }, ({ event }) => {
        this.#live(run, () => {
          this.#backend.appendRunEvent(run.runId, event);
        });
      })
      .with({ kind: "output-begin" }, ({ path }) => {
        run.state = "outputs";
        this.#live(run, () => {
          this.#backend.beginRunOutput(run.runId, path);
        });
      })
      .with({ kind: "output-chunk" }, ({ path, data }) => {
        this.#live(run, () => {
          this.#backend.appendRunOutput(run.runId, path, data);
        });

        if (this.#current(run, attempt)) {
          this.#post(run, { kind: "output-ack", attempt });
        }
      })
      .with({ kind: "done" }, ({ outcome }) => {
        this.#complete(run, outcome);
      })
      .exhaustive();
  }

  async #sendInputs(run: RunJob, attempt: number): Promise<void> {
    for (const path of run.inputPaths) {
      if (!(await this.#sendInput(run, attempt, path))) {
        return;
      }
    }

    if (this.#current(run, attempt)) {
      run.state = "computing";
      this.#post(run, { kind: "run", attempt });
    }
  }

  async #sendInput(run: RunJob, attempt: number, path: string): Promise<boolean> {
    let reader: WasmChunkReader;

    try {
      reader = this.#backend.readRunInput(run.runId, path);
    } catch (error) {
      this.#liveFailed(run, error);

      return false;
    }

    try {
      this.#post(run, { kind: "input-begin", attempt, path });

      for (;;) {
        if (!(await run.window.acquire()) || !this.#current(run, attempt)) {
          return false;
        }

        const chunk = reader.chunk();

        if (chunk === undefined) {
          run.window.release();

          return true;
        }

        this.#post(run, { kind: "input-chunk", attempt, path, data: chunk }, transferList(chunk));
      }
    } catch (error) {
      this.#liveFailed(run, error);

      return false;
    } finally {
      reader.free();
    }
  }

  #lost(run: RunJob, attempt: number, lines: string[]): void {
    if (!this.#current(run, attempt)) {
      return;
    }

    if (run.state === "starting") {
      this.#startFailed(run, attempt, lines);
    } else {
      this.#fail(run, lines);
    }
  }

  #startFailed(run: RunJob, attempt: number, lines: string[]): void {
    if (!this.#current(run, attempt) || run.state !== "starting") {
      return;
    }

    if (run.variant === "threaded") {
      this.#threadedFailed = true;
      this.#stop(run);
      console.warn(`[TreeTime] the threaded run worker could not start; runs use one thread: ${lines.join("\n")}`);
      this.#launch(run);
    } else {
      this.#fail(run, lines);
    }
  }

  #cancel(runId: string): void {
    const run = this.#jobs.get(runId);

    if (run !== undefined) {
      this.#fail(run, []);
    }
  }

  #complete(run: RunJob, outcome: string): void {
    this.#settle(run);

    try {
      this.#backend.completeRun(run.runId, outcome);
    } catch (error) {
      if (!this.#trap(error)) {
        this.#report(run, [errorMessage(error)]);
      }
    }
  }

  #fail(run: RunJob, lines: string[]): void {
    this.#settle(run);
    this.#report(run, lines);
  }

  #report(run: RunJob, lines: string[]): void {
    try {
      this.#backend.failRun(run.runId, lines);
    } catch (error) {
      if (!this.#trap(error)) {
        console.error(`[TreeTime] the end of run ${run.runId} could not be recorded`, error);
      }
    }
  }

  #live(run: RunJob, work: () => void): void {
    try {
      work();
    } catch (error) {
      this.#liveFailed(run, error);
    }
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- a thrown value of the back end can be anything; its message is read from it
  #liveFailed(run: RunJob, error: unknown): void {
    if (run.state !== "settled" && !this.#trap(error)) {
      this.#fail(run, [errorMessage(error)]);
    }
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- a thrown value of the back end can be anything
  #backendFailed(error: unknown): void {
    if (!this.#trap(error)) {
      console.error("[TreeTime] the run orders of the back end could not be read", error);
    }
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- a thrown value of the back end can be anything
  #trap(error: unknown): boolean {
    if (!(error instanceof WebAssembly.RuntimeError)) {
      return false;
    }

    if (!this.#trapped) {
      this.#trapped = true;
      [...this.#jobs.values()].forEach((run) => {
        this.#settle(run);
      });
      this.#onTrap(error);
    }

    return true;
  }

  #settle(run: RunJob): void {
    this.#stop(run);
    run.state = "settled";
    this.#jobs.delete(run.runId);
  }

  #stop(run: RunJob): void {
    clearTimeout(run.timer);
    run.window.close();
    run.worker?.terminate();
    run.worker = undefined;
  }

  #current(run: RunJob, attempt: number): boolean {
    return run.attempt === attempt && run.state !== "settled";
  }

  #post(run: RunJob, message: HostMessage, transfer: Transferable[] = []): void {
    run.worker?.postMessage(message, transfer);
  }
}

export interface JobBackend {
  takeJobOrders(): WasmJobOrder[];
  readRunInput(runId: string, path: string): WasmChunkReader;
  appendRunEvent(runId: string, eventJson: string): void;
  beginRunOutput(runId: string, relativePath: string): void;
  appendRunOutput(runId: string, relativePath: string, chunk: Uint8Array): void;
  completeRun(runId: string, outcomeJson: string): boolean;
  failRun(runId: string, lines: string[]): boolean;
}

export interface WasmJobOrder {
  readonly kind: string;
  readonly runId: string;
  readonly job?: string | undefined;
  readonly inputPaths: string[];
  free(): void;
}

export interface JobModules {
  single: WebAssembly.Module;
  threaded?: ThreadedModule | undefined;
}

export interface ThreadedModule {
  module(): Promise<WebAssembly.Module>;
  jobs: number;
}

export interface JobWorker {
  postMessage(message: HostMessage, transfer: Transferable[]): void;
  addEventListener(type: "message", listener: (event: { data: unknown }) => void): void;
  addEventListener(type: "error", listener: (event: { message?: string | undefined }) => void): void;
  terminate(): void;
}

function readOrder(order: WasmJobOrder): JobOrder {
  const { kind, runId, job, inputPaths } = order;

  order.free();

  if (kind === "cancel") {
    return { kind, runId };
  }

  if (kind !== "start" || job === undefined) {
    throw new Error(`the back end queued an unknown run order \`${kind}\` for run \`${runId}\``);
  }

  return { kind, runId, job, inputPaths };
}

type JobOrder = { kind: "start"; runId: string; job: string; inputPaths: string[] } | { kind: "cancel"; runId: string };

interface RunJob {
  readonly runId: string;
  readonly job: string;
  readonly inputPaths: string[];
  state: JobState;
  attempt: number;
  variant: ModuleVariant;
  worker: JobWorker | undefined;
  window: ChunkWindow;
  timer: ReturnType<typeof setTimeout> | undefined;
}

type JobState = "starting" | "inputs" | "computing" | "outputs" | "settled";
