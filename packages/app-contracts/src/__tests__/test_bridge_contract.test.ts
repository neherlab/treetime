import { describe, expect, test } from "vitest";
import { ZodError } from "zod";

import {
  CancelledError,
  CommandError,
  RunEndedError,
  createBridge,
  parseRunEvent,
  type BridgeTransport,
  type CommandOutcome,
  type LogEvent,
  type ProgressEvent,
  type TransportEventOptions,
} from "../index";

type ParsedRunEvent = ReturnType<typeof parseRunEvent>;

const OUTCOME: CommandOutcome = {
  command: "ancestral",
  output_files: [{ path: "/runs/r1/out/ancestral.nwk", kind: "nwk" }],
};

const RECORD = {
  id: "r1",
  title: "ancestral",
  command: "ancestral",
  config: { tree: "t.nwk" },
  status: "running",
  pinned: false,
  created_at: "2026-09-25T10:00:00Z",
  started_at: "2026-09-25T10:00:01Z",
  finished_at: null,
  duration_seconds: null,
  treetime_version: "1.0.0",
  inputs: [],
  config_hash: null,
  changed_settings: [],
  headline: {},
  output_files: [],
  error: null,
};

function stubTransport(overrides: Partial<BridgeTransport>): BridgeTransport {
  const missing = (name: string) => () => Promise.reject(new Error(`${name} not stubbed`));

  return {
    version: missing("version"),
    datasets: missing("datasets"),
    checkConfig: missing("checkConfig"),
    runConfig: missing("runConfig"),
    checkInputs: missing("checkInputs"),
    listRuns: missing("listRuns"),
    createRun: missing("createRun"),
    getRun: missing("getRun"),
    startRun: missing("startRun"),
    updateRun: missing("updateRun"),
    cancelRun: missing("cancelRun"),
    deleteRun: missing("deleteRun"),
    restoreRun: missing("restoreRun"),
    purgeRun: missing("purgeRun"),
    runEvents: missing("runEvents"),
    runFiles: missing("runFiles"),
    readRunFile: missing("readRunFile"),
    runArchive: missing("runArchive"),
    uploadInput: missing("uploadInput"),
    ...overrides,
  };
}

function events(...items: Array<{ type: string; data: unknown }>): ParsedRunEvent[] {
  return items.map((item, seq) =>
    parseRunEvent({ seq, time: "2026-09-25T10:00:00Z", type: item.type, data: item.data }),
  );
}

function replay(list: ParsedRunEvent[]): BridgeTransport["runEvents"] {
  return (_id: string, options: TransportEventOptions) => {
    for (const event of list.filter((item) => item.seq >= options.from)) {
      options.onEvent(event);
    }

    return Promise.resolve();
  };
}

describe("bridge result validation", () => {
  test("version returns the validated result", async () => {
    const bridge = createBridge(stubTransport({ version: () => Promise.resolve({ version: "1.2.3" }) }));
    await expect(bridge.version()).resolves.toStrictEqual({ version: "1.2.3" });
  });

  test("datasets returns the datasets and the example configurations", async () => {
    const catalog = {
      data_dir: "data",
      datasets: [{ name: "zika/20", files: ["tree.nwk"] }],
      examples: [{ path: "zika/20/mugration.yaml", command: "mugration", title: "Geography", content: "tree: x" }],
    };

    const bridge = createBridge(stubTransport({ datasets: () => Promise.resolve(catalog) }));
    await expect(bridge.datasets()).resolves.toStrictEqual(catalog);
  });

  test("a malformed result rejects with a ZodError", async () => {
    const bridge = createBridge(stubTransport({ version: () => Promise.resolve({ version: 1 }) }));
    await expect(bridge.version()).rejects.toBeInstanceOf(ZodError);
  });

  test("checkConfig passes the request through and validates the response", async () => {
    let captured: unknown;
    const response = { status: "valid", config: { tree: "t.nwk", max_iter: 2 } };

    const bridge = createBridge(
      stubTransport({
        checkConfig: (request) => {
          captured = request;

          return Promise.resolve(response);
        },
      }),
    );

    await expect(bridge.checkConfig({ command: "timetree", text: "tree: t.nwk\n" })).resolves.toStrictEqual(response);
    expect(captured).toStrictEqual({ command: "timetree", text: "tree: t.nwk\n" });
  });

  test("runConfig passes the request through and validates the response", async () => {
    let captured: unknown;

    const response = {
      status: "valid",
      config: { tree: "t.nwk", output_all: "out", output_selection: ["Auspice"] },
      config_hash: "abc",
      config_hash_error: null,
    };

    const bridge = createBridge(
      stubTransport({
        runConfig: (request) => {
          captured = request;

          return Promise.resolve(response);
        },
      }),
    );

    await expect(bridge.runConfig({ command: "prune", config: { tree: "t.nwk" } })).resolves.toStrictEqual(response);
    expect(captured).toStrictEqual({ command: "prune", config: { tree: "t.nwk" } });
  });

  test("checkInputs validates the facts", async () => {
    const facts = {
      tree: { tips: 3, internal_nodes: 1, polytomies: 1, unnamed_tips: 0, duplicate_tip_names: [] },
      alignment: null,
      metadata: null,
      tips_without_metadata: null,
      tips_without_sequence: null,
      problems: [{ input: "metadata", message: "cannot read" }],
    };

    const bridge = createBridge(stubTransport({ checkInputs: () => Promise.resolve(facts) }));
    await expect(bridge.checkInputs({ tree: "t.nwk" })).resolves.toStrictEqual(facts);
  });

  test("listRuns validates the run list", async () => {
    const list = { runs: [], active_runs: 2 };
    const bridge = createBridge(stubTransport({ listRuns: () => Promise.resolve(list) }));
    await expect(bridge.listRuns()).resolves.toStrictEqual(list);
  });

  test("cancelRun returns whether cancellation was requested", async () => {
    const bridge = createBridge(stubTransport({ cancelRun: () => Promise.resolve({ cancelled: true }) }));
    await expect(bridge.cancelRun("r1")).resolves.toBe(true);
  });
});

describe("bridge run events", () => {
  test("followRun resolves with the terminal event and forwards every event", async () => {
    const stream = events(
      { type: "started", data: { job_id: "r1", command: "clock" } },
      { type: "progress", data: { stage: "read", fraction: 0.5, message: "" } },
      { type: "terminal", data: { status: "cancelled", job_id: "r1" } },
    );

    const seen: number[] = [];
    const bridge = createBridge(stubTransport({ runEvents: replay(stream) }));

    const terminal = await bridge.followRun("r1", {
      onEvent: (event) => {
        seen.push(event.seq);
      },
    });

    expect(terminal).toStrictEqual({ status: "cancelled", job_id: "r1" });
    expect(seen).toStrictEqual([0, 1, 2]);
  });

  test("followRun resumes from an offset", async () => {
    const stream = events(
      { type: "started", data: { job_id: "r1", command: "clock" } },
      { type: "log", data: { level: "info", message: "a" } },
      { type: "terminal", data: { status: "cancelled", job_id: "r1" } },
    );

    const seen: number[] = [];
    const bridge = createBridge(stubTransport({ runEvents: replay(stream) }));
    await bridge.followRun("r1", {
      from: 1,
      onEvent: (event) => {
        seen.push(event.seq);
      },
    });
    expect(seen).toStrictEqual([1, 2]);
  });

  test("followRun rejects when the stream ends without a terminal event", async () => {
    const stream = events({ type: "started", data: { job_id: "r1", command: "clock" } });
    const bridge = createBridge(stubTransport({ runEvents: replay(stream) }));
    await expect(bridge.followRun("r1")).rejects.toBeInstanceOf(RunEndedError);
  });

  test("parseRunEvent reads non-finite iteration values", () => {
    const event = {
      seq: 7,
      time: "2026-09-25T10:00:00Z",
      type: "iteration",
      data: {
        iteration: 2,
        n_diff: 0,
        n_resolved: 0,
        max_time_change: 0.1,
        rms_time_change: 0.01,
        log_lh_seq: null,
        log_lh_pos: -12.5,
        log_lh_coal: "inf",
        log_lh_total: "inf",
        clock_rate: 0.001,
        r_squared: 0.8,
      },
    };

    expect(parseRunEvent(event)).toStrictEqual(event);
  });

  test("parseRunEvent rejects an unknown event type", () => {
    expect(() => parseRunEvent({ seq: 0, time: "t", type: "result", data: {} })).toThrow(ZodError);
  });

  test("parseRunEvent rejects an event without its position", () => {
    expect(() => parseRunEvent({ type: "log", data: { level: "info", message: "x" } })).toThrow(ZodError);
  });
});

describe("bridge commands", () => {
  test("a command creates a run, follows its events and resolves with the outcome", async () => {
    const stream = events(
      { type: "started", data: { job_id: "r1", command: "ancestral" } },
      { type: "progress", data: { stage: "read", fraction: 0.25, message: "reading" } },
      { type: "log", data: { level: "warn", message: "tip without date" } },
      {
        type: "iteration",
        data: { iteration: 1, n_diff: 3, n_resolved: 0, clock_rate: 0.001 },
      },
      { type: "terminal", data: { status: "ok", job_id: "r1", result: OUTCOME } },
    );

    let created: unknown;
    const started: string[] = [];
    const progress: ProgressEvent[] = [];
    const logs: LogEvent[] = [];
    const iterations: number[] = [];

    const bridge = createBridge(
      stubTransport({
        createRun: (request) => {
          created = request;

          return Promise.resolve(RECORD);
        },
        runEvents: replay(stream),
      }),
    );

    const outcome = await bridge.ancestral(
      { tree: "t.nwk" },
      {
        title: "zika",
        onStarted: (id) => {
          started.push(id);
        },
        onProgress: (event) => {
          progress.push(event);
        },
        onLog: (event) => {
          logs.push(event);
        },
        onIteration: (event) => {
          iterations.push(event.iteration);
        },
      },
    );

    expect(outcome).toStrictEqual(OUTCOME);
    expect(created).toStrictEqual({
      command: "ancestral",
      config: { tree: "t.nwk" },
      defer_start: false,
      title: "zika",
    });
    expect(started).toStrictEqual(["r1"]);
    expect(progress).toStrictEqual([{ stage: "read", fraction: 0.25, message: "reading" }]);
    expect(logs).toStrictEqual([{ level: "warn", message: "tip without date" }]);
    expect(iterations).toStrictEqual([1]);
  });

  test("an error terminal event rejects with a CommandError carrying the cause chain", async () => {
    const stream = events({
      type: "terminal",
      data: { status: "error", job_id: "r1", message: "When reading dates", causes: ["file not found"] },
    });

    const bridge = createBridge(stubTransport({ createRun: () => Promise.resolve(RECORD), runEvents: replay(stream) }));
    const error: unknown = await bridge.clock({ tree: "t.nwk" }).catch((err: unknown) => err);

    expect(error).toBeInstanceOf(CommandError);
    expect(error).toMatchObject({ jobId: "r1", message: "When reading dates", causes: ["file not found"] });
  });

  test("a cancelled terminal event rejects with a CancelledError", async () => {
    const stream = events({ type: "terminal", data: { status: "cancelled", job_id: "r1" } });
    const bridge = createBridge(stubTransport({ createRun: () => Promise.resolve(RECORD), runEvents: replay(stream) }));
    await expect(bridge.timetree({})).rejects.toBeInstanceOf(CancelledError);
  });

  test("an interrupted terminal event rejects with a CommandError", async () => {
    const stream = events({ type: "terminal", data: { status: "interrupted", job_id: "r1" } });
    const bridge = createBridge(stubTransport({ createRun: () => Promise.resolve(RECORD), runEvents: replay(stream) }));
    await expect(bridge.prune({ tree: "t.nwk" })).rejects.toBeInstanceOf(CommandError);
  });

  test("aborting the signal requests cancellation of the run", async () => {
    const controller = new AbortController();
    const cancelled: string[] = [];

    const bridge = createBridge(
      stubTransport({
        createRun: () => Promise.resolve(RECORD),
        cancelRun: (id) => {
          cancelled.push(id);

          return Promise.resolve({ cancelled: true });
        },
        runEvents: (_id, options) => {
          options.onEvent({ seq: 0, time: "t", type: "started", data: { job_id: "r1", command: "optimize" } });
          controller.abort();
          options.onEvent({ seq: 1, time: "t", type: "terminal", data: { status: "cancelled", job_id: "r1" } });

          return Promise.resolve();
        },
      }),
    );

    await expect(bridge.optimize({ tree: "t.nwk" }, { signal: controller.signal })).rejects.toBeInstanceOf(
      CancelledError,
    );
    expect(cancelled).toStrictEqual(["r1"]);
  });

  test("an already aborted signal rejects before creating a run", async () => {
    const controller = new AbortController();
    controller.abort();
    const bridge = createBridge(stubTransport({}));
    await expect(bridge.clock({}, { signal: controller.signal })).rejects.toBeInstanceOf(CancelledError);
  });

  test("cancelled error has a stable Error name", () => {
    const err = new CancelledError();
    expect(err).toBeInstanceOf(Error);
    expect(err.name).toBe("CancelledError");
  });
});
