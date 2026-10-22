import pWaitFor from "p-wait-for";

import type { HostMessage } from "../../job/protocol";
import type { WasmChunkReader } from "../../job/serve";
import { type JobBackend, JobHost, type JobWorker, type ThreadedModule, type WasmJobOrder } from "../jobs";

const SINGLE_MODULE = new WebAssembly.Module(new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]));

const THREADED_MODULE = new WebAssembly.Module(
  new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00]),
);

const JOB = '{"id":"r1"}';

const TREE = "/treetime/runs/r1/inputs/1/tree.nwk";

const ALIGNMENT = "/treetime/runs/r1/inputs/2/aln.fasta";

const SIX_CHUNKS = { [TREE]: ["1", "2", "3", "4", "5", "6"] };

const BLOCKED_WORKER = "SecurityError: Failed to construct 'Worker': the script is blocked";

describe("job host start", () => {
  test("a start order starts a run worker with the single-threaded module, the job, and one thread", async () => {
    const { host, workers } = fakeHost();

    host.backend.orders.push(startOrder());
    host.jobs.drain();
    await posted(workers, 0, 1);

    expect([workers.length, workers[0]?.posted]).toStrictEqual([
      1,
      [{ kind: "start", attempt: 1, variant: "single", module: SINGLE_MODULE, job: JOB, jobs: 1 }],
    ]);
  });

  test("a started worker receives every input file in chunks and then the run message", async () => {
    const { host, workers } = fakeHost({ inputs: { [TREE]: ["(A,", "B);"], [ALIGNMENT]: [] } });

    await startedJob(host, workers);
    await posted(workers, 0, 6);

    expect(workers[0]?.posted.slice(1).map(describeMessage)).toStrictEqual([
      `input-begin ${TREE}`,
      `input-chunk ${TREE} (A,`,
      `input-chunk ${TREE} B);`,
      `input-begin ${ALIGNMENT}`,
      "run",
    ]);
  });

  test("the input chunks travel as transferred buffers", async () => {
    const { host, workers } = fakeHost({ inputs: { [TREE]: ["(A,B);"] } });

    await startedJob(host, workers, [TREE]);
    await posted(workers, 0, 4);

    const transferred = workers[0]?.transferred[2]?.map((buffer) => buffer instanceof ArrayBuffer && buffer.byteLength);

    expect([describeMessage(workers[0]?.posted[2]), transferred]).toStrictEqual([`input-chunk ${TREE} (A,B);`, [6]]);
  });

  test("a receiver that acknowledges slowly holds at most 4 chunks without an acknowledgement", async () => {
    const { host, workers } = fakeHost({ inputs: SIX_CHUNKS });

    await startedJob(host, workers, [TREE]);
    await posted(workers, 0, 6);
    const beforeAck = countChunks(workers[0]);
    workers[0]?.emit({ kind: "input-ack", attempt: 1 });
    await posted(workers, 0, 7);
    const afterAck = countChunks(workers[0]);

    expect([beforeAck, afterAck]).toStrictEqual([4, 5]);
  });
});

describe("job host run", () => {
  test("events reach the back end in the order the worker sends them", async () => {
    const { host, workers } = fakeHost();

    await computingJob(host, workers);
    workers[0]?.emit({ kind: "event", attempt: 1, event: '{"type":"log","n":1}' });
    workers[0]?.emit({ kind: "event", attempt: 1, event: '{"type":"log","n":2}' });

    expect(host.backend.calls).toStrictEqual([
      ["appendRunEvent", "r1", '{"type":"log","n":1}'],
      ["appendRunEvent", "r1", '{"type":"log","n":2}'],
    ]);
  });

  test("output files are written, each chunk is acknowledged, and done completes the run and ends the worker", async () => {
    const { host, workers } = fakeHost();
    const chunk = new TextEncoder().encode("(A,B);");

    await computingJob(host, workers);
    workers[0]?.emit({ kind: "output-begin", attempt: 1, path: "empty.txt" });
    workers[0]?.emit({ kind: "output-begin", attempt: 1, path: "ancestral.nwk" });
    workers[0]?.emit({ kind: "output-chunk", attempt: 1, path: "ancestral.nwk", data: chunk });
    workers[0]?.emit({ kind: "done", attempt: 1, outcome: '{"terminal":{}}' });

    expect([host.backend.calls, workers[0]?.posted.at(-1), workers[0]?.terminated]).toStrictEqual([
      [
        ["beginRunOutput", "r1", "empty.txt"],
        ["beginRunOutput", "r1", "ancestral.nwk"],
        ["appendRunOutput", "r1", "ancestral.nwk", chunk],
        ["completeRun", "r1", '{"terminal":{}}'],
      ],
      { kind: "output-ack", attempt: 1 },
      true,
    ]);
  });

  test("done followed by a cancel completes the run once and records no failure", async () => {
    const { host, workers } = fakeHost();

    await computingJob(host, workers);
    workers[0]?.emit({ kind: "done", attempt: 1, outcome: "{}" });
    host.backend.orders.push(cancelOrder());
    host.jobs.drain();

    expect(host.backend.calls).toStrictEqual([["completeRun", "r1", "{}"]]);
  });

  test("a crash ends the run with the error lines of the worker", async () => {
    const { host, workers } = fakeHost();

    await computingJob(host, workers);
    workers[0]?.emit({ kind: "crashed", attempt: 1, lines: ["memory allocation of 8 bytes failed", "Aborted()"] });

    expect([host.backend.calls, workers[0]?.terminated]).toStrictEqual([
      [["failRun", "r1", ["memory allocation of 8 bytes failed", "Aborted()"]]],
      true,
    ]);
  });

  test("a crash and an error event of the same abort end the run once", async () => {
    const { host, workers } = fakeHost();

    await computingJob(host, workers);
    workers[0]?.emit({ kind: "crashed", attempt: 1, lines: ["stack overflow"] });
    workers[0]?.fail("RuntimeError: Aborted(stack overflow)");

    expect(host.backend.calls).toStrictEqual([["failRun", "r1", ["stack overflow"]]]);
  });

  test("an error event of a running worker ends the run with its message", async () => {
    const { host, workers } = fakeHost();

    await computingJob(host, workers);
    workers[0]?.fail("Uncaught ZodError");

    expect(host.backend.calls).toStrictEqual([["failRun", "r1", ["Uncaught ZodError"]]]);
  });

  test("a back end that rejects an output path ends the run with its error", async () => {
    const { host, workers } = fakeHost({
      rejectOutput: "file path `../run.json` must name a file inside the run's output folder",
    });

    await computingJob(host, workers);
    workers[0]?.emit({ kind: "output-begin", attempt: 1, path: "../run.json" });

    expect([host.backend.calls.at(-1), workers[0]?.terminated]).toStrictEqual([
      ["failRun", "r1", ["file path `../run.json` must name a file inside the run's output folder"]],
      true,
    ]);
  });
});

describe("job host cancel", () => {
  test.each([
    ["starting", startingJob, undefined],
    ["inputs", inputsJob, SIX_CHUNKS],
    ["computing", computingJob, undefined],
    ["outputs", outputsJob, undefined],
  ])("a cancel while the job is %s ends the worker and records the end", async (_state, reach, inputs) => {
    const { host, workers } = fakeHost(inputs === undefined ? {} : { inputs });

    await reach(host, workers);
    const calls = host.backend.calls.length;
    host.backend.orders.push(cancelOrder());
    host.jobs.drain();

    expect([host.backend.calls.slice(calls), workers[0]?.terminated]).toStrictEqual([[["failRun", "r1", []]], true]);
  });

  test("a cancel stops an input transfer between chunks", async () => {
    const { host, workers } = fakeHost({ inputs: SIX_CHUNKS });

    await inputsJob(host, workers);
    host.backend.orders.push(cancelOrder());
    host.jobs.drain();
    workers[0]?.emit({ kind: "input-ack", attempt: 1 });
    await pWaitFor(() => host.backend.readersFreed === 1);

    expect([countChunks(workers[0]), host.backend.readersFreed]).toStrictEqual([4, 1]);
  });
});

describe("job host start failures", () => {
  test("a start failure of the threaded module starts the job again on the single-threaded module, also for later jobs", async () => {
    const { host, workers } = fakeHost({ threaded: { module: () => Promise.resolve(THREADED_MODULE), jobs: 8 } });

    host.backend.orders.push(startOrder());
    host.jobs.drain();
    await posted(workers, 0, 1);
    workers[0]?.emit({ kind: "failed", attempt: 1, lines: ["Aborted(pthread worker failed to load)"] });
    await posted(workers, 1, 1);
    host.backend.orders.push(startOrder("r2"));
    host.jobs.drain();
    await posted(workers, 2, 1);

    expect([workers.map(({ posted: [start] }) => start), workers[0]?.terminated, host.backend.calls]).toStrictEqual([
      [
        { kind: "start", attempt: 1, variant: "threaded", module: THREADED_MODULE, job: JOB, jobs: 8 },
        { kind: "start", attempt: 2, variant: "single", module: SINGLE_MODULE, job: JOB, jobs: 1 },
        { kind: "start", attempt: 3, variant: "single", module: SINGLE_MODULE, job: JOB, jobs: 1 },
      ],
      true,
      [],
    ]);
  });

  test("a threaded module that does not compile is a start failure of the threaded attempt", async () => {
    const { host, workers } = fakeHost({
      threaded: { module: () => Promise.reject(new Error("HTTP status 404")), jobs: 8 },
    });

    host.backend.orders.push(startOrder());
    host.jobs.drain();
    await posted(workers, 1, 1);

    expect([workers[0]?.posted, workers[1]?.posted[0]?.kind]).toStrictEqual([[], "start"]);
  });

  test("messages of a replaced attempt are ignored", async () => {
    const { host, workers } = fakeHost({ threaded: { module: () => Promise.resolve(THREADED_MODULE), jobs: 8 } });

    host.backend.orders.push(startOrder());
    host.jobs.drain();
    await posted(workers, 0, 1);
    workers[0]?.fail("Uncaught (in promise) Error: worker failed to load");
    await posted(workers, 1, 1);
    workers[0]?.emit({ kind: "started", attempt: 1 });
    workers[0]?.emit({ kind: "crashed", attempt: 1, lines: ["late"] });
    workers[1]?.emit({ kind: "crashed", attempt: 1, lines: ["of attempt 1"] });

    expect([workers[1]?.posted.length, workers[1]?.terminated, host.backend.calls]).toStrictEqual([1, false, []]);
  });

  test("a start failure of the single-threaded module ends the run with the error lines", async () => {
    const { host, workers } = fakeHost();

    host.backend.orders.push(startOrder());
    host.jobs.drain();
    await posted(workers, 0, 1);
    workers[0]?.emit({ kind: "failed", attempt: 1, lines: ["CompileError: wasm validation error"] });

    expect([host.backend.calls, workers.length]).toStrictEqual([
      [["failRun", "r1", ["CompileError: wasm validation error"]]],
      1,
    ]);
  });

  test("a threaded worker that cannot be constructed starts the job again on the single-threaded module", async () => {
    const { host, workers } = fakeHost({
      threaded: { module: () => Promise.resolve(THREADED_MODULE), jobs: 8 },
      throwingStarts: 1,
    });

    host.backend.orders.push(startOrder());
    host.jobs.drain();
    await posted(workers, 0, 1);

    expect([workers.map(({ posted: [start] }) => start), host.backend.calls]).toStrictEqual([
      [{ kind: "start", attempt: 2, variant: "single", module: SINGLE_MODULE, job: JOB, jobs: 1 }],
      [],
    ]);
  });

  test("a single-threaded worker that cannot be constructed ends the run with the error", () => {
    const { host, workers } = fakeHost({ throwingStarts: 1 });

    host.backend.orders.push(startOrder());
    host.jobs.drain();

    expect([host.backend.calls, workers.length]).toStrictEqual([[["failRun", "r1", [BLOCKED_WORKER]]], 0]);
  });

  test("a worker that does not start in time ends the run", async () => {
    const { host, workers } = fakeHost({ startTimeoutMs: 1 });

    host.backend.orders.push(startOrder());
    host.jobs.drain();
    await pWaitFor(() => host.backend.calls.length === 1);

    expect([host.backend.calls, workers[0]?.terminated]).toStrictEqual([
      [["failRun", "r1", ["the run worker did not start within 0.001 s"]]],
      true,
    ]);
  });
});

describe("job host back end failures", () => {
  test("a trap of the back end ends every job and is reported once", async () => {
    const traps: string[] = [];

    const { host, workers } = fakeHost({
      trapOnEvent: true,
      onTrap: (trap) => {
        traps.push(trap.message);
      },
    });

    await computingJob(host, workers);
    workers[0]?.emit({ kind: "event", attempt: 1, event: "{}" });
    host.jobs.drain();

    expect([traps, workers[0]?.terminated, host.backend.calls]).toStrictEqual([
      ["unreachable"],
      true,
      [["appendRunEvent", "r1", "{}"]],
    ]);
  });
});

interface FakeSetup {
  host: { jobs: JobHost; backend: FakeBackend };
  workers: FakeWorker[];
}

interface FakeBackend extends JobBackend {
  orders: WasmJobOrder[];
  calls: unknown[][];
  readonly readersFreed: number;
}

interface FakeWorker extends JobWorker {
  posted: HostMessage[];
  transferred: Transferable[][];
  terminated: boolean;
  emit(data: unknown): void;
  fail(message: string): void;
}

function fakeHost({
  inputs = { [TREE]: ["(A,B);"] },
  threaded,
  startTimeoutMs,
  rejectOutput,
  throwingStarts = 0,
  trapOnEvent = false,
  onTrap = () => undefined,
}: {
  inputs?: Record<string, string[]>;
  threaded?: ThreadedModule;
  startTimeoutMs?: number;
  rejectOutput?: string;
  throwingStarts?: number;
  trapOnEvent?: boolean;
  onTrap?: (trap: WebAssembly.RuntimeError) => void;
} = {}): FakeSetup {
  const workers: FakeWorker[] = [];
  const backend = fakeBackend(inputs, { rejectOutput, trapOnEvent });
  let starts = 0;

  const jobs = new JobHost({
    backend,
    startWorker: () => {
      starts += 1;

      if (starts <= throwingStarts) {
        throw new Error(BLOCKED_WORKER);
      }

      const worker = fakeWorker();

      workers.push(worker);

      return worker;
    },
    modules: { single: SINGLE_MODULE, threaded },
    onTrap,
    startTimeoutMs,
  });

  return { host: { jobs, backend }, workers };
}

function fakeBackend(
  inputs: Record<string, string[]>,
  { rejectOutput, trapOnEvent }: { rejectOutput: string | undefined; trapOnEvent: boolean },
): FakeBackend {
  let readersFreed = 0;

  const backend: FakeBackend = {
    orders: [],
    calls: [],
    get readersFreed() {
      return readersFreed;
    },
    takeJobOrders: () => backend.orders.splice(0),
    readRunInput: (_runId, path): WasmChunkReader => {
      const chunks = (inputs[path] ?? []).map((text) => new TextEncoder().encode(text));

      return {
        chunk: () => chunks.shift(),
        free: () => {
          readersFreed += 1;
        },
      };
    },
    appendRunEvent: (runId, eventJson) => {
      backend.calls.push(["appendRunEvent", runId, eventJson]);

      if (trapOnEvent) {
        throw new WebAssembly.RuntimeError("unreachable");
      }
    },
    beginRunOutput: (runId, relativePath) => {
      backend.calls.push(["beginRunOutput", runId, relativePath]);

      if (rejectOutput !== undefined) {
        throw new Error(rejectOutput);
      }
    },
    appendRunOutput: (runId, relativePath, chunk) => {
      backend.calls.push(["appendRunOutput", runId, relativePath, chunk]);
    },
    completeRun: (runId, outcomeJson) => {
      backend.calls.push(["completeRun", runId, outcomeJson]);

      return true;
    },
    failRun: (runId, lines) => {
      backend.calls.push(["failRun", runId, lines]);

      return true;
    },
  };

  return backend;
}

function fakeWorker(): FakeWorker {
  const messageListeners: Array<(event: { data: unknown }) => void> = [];
  const errorListeners: Array<(event: { message?: string | undefined }) => void> = [];

  const worker: FakeWorker = {
    posted: [],
    transferred: [],
    terminated: false,
    postMessage(message, transfer) {
      worker.posted.push(message);
      worker.transferred.push(transfer);
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

async function startingJob(host: FakeSetup["host"], workers: FakeWorker[]): Promise<void> {
  host.backend.orders.push(startOrder());
  host.jobs.drain();
  await posted(workers, 0, 1);
}

async function inputsJob(host: FakeSetup["host"], workers: FakeWorker[]): Promise<void> {
  await startedJob(host, workers, [TREE]);
  await posted(workers, 0, 6);
}

async function outputsJob(host: FakeSetup["host"], workers: FakeWorker[]): Promise<void> {
  await computingJob(host, workers);
  workers[0]?.emit({ kind: "output-begin", attempt: 1, path: "ancestral.nwk" });
}

async function computingJob(host: FakeSetup["host"], workers: FakeWorker[]): Promise<void> {
  await startedJob(host, workers, [TREE]);
  await pWaitFor(() => workers[0]?.posted.at(-1)?.kind === "run");
}

async function startedJob(
  host: FakeSetup["host"],
  workers: FakeWorker[],
  inputPaths = [TREE, ALIGNMENT],
): Promise<void> {
  host.backend.orders.push({ ...startOrder(), inputPaths });
  host.jobs.drain();
  await posted(workers, 0, 1);
  workers[0]?.emit({ kind: "started", attempt: 1 });
}

function startOrder(runId = "r1"): WasmJobOrder {
  return { kind: "start", runId, job: JOB, inputPaths: [TREE], free: () => undefined };
}

function cancelOrder(): WasmJobOrder {
  return { kind: "cancel", runId: "r1", job: undefined, inputPaths: [], free: () => undefined };
}

async function posted(workers: FakeWorker[], index: number, count: number): Promise<void> {
  await pWaitFor(() => (workers[index]?.posted.length ?? 0) >= count);
}

function countChunks(worker: FakeWorker | undefined): number {
  return worker?.posted.filter(({ kind }) => kind === "input-chunk").length ?? 0;
}

function describeMessage(message: HostMessage | undefined): string | undefined {
  if (message?.kind === "input-chunk") {
    return `input-chunk ${message.path} ${new TextDecoder().decode(message.data)}`;
  }

  return message?.kind === "input-begin" ? `input-begin ${message.path}` : message?.kind;
}
