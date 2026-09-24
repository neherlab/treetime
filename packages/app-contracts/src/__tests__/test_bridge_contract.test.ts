import { describe, expect, test } from "vitest";
import { ZodError } from "zod";

import {
  CancelledError,
  CommandError,
  createBridge,
  parseJobEvent,
  type BridgeTransport,
  type JobEvent,
  type LogEvent,
  type ProgressEvent,
  type TransportCommandOptions,
} from "../index";

function stubTransport(overrides: Partial<BridgeTransport>): BridgeTransport {
  return {
    query: overrides.query ?? (() => Promise.reject(new Error("query not stubbed"))),
    request: overrides.request ?? (() => Promise.reject(new Error("request not stubbed"))),
    command: overrides.command ?? (() => Promise.reject(new Error("command not stubbed"))),
  };
}

const OUTCOME = { command: "ancestral", output_files: ["out/ancestral.nwk"] };

describe("bridge result validation", () => {
  test("version returns the validated result", async () => {
    const bridge = createBridge(stubTransport({ query: () => Promise.resolve({ version: "1.2.3" }) }));
    await expect(bridge.version()).resolves.toStrictEqual({ version: "1.2.3" });
  });

  test("datasets validates and returns an array", async () => {
    const datasets = [{ name: "flu", files: ["tree.nwk"] }];
    const bridge = createBridge(stubTransport({ query: () => Promise.resolve(datasets) }));
    await expect(bridge.datasets()).resolves.toStrictEqual(datasets);
  });

  test("a malformed query result rejects with a ZodError", async () => {
    const bridge = createBridge(stubTransport({ query: () => Promise.resolve({ version: 1 }) }));
    await expect(bridge.version()).rejects.toBeInstanceOf(ZodError);
  });

  test("checkConfig sends the request to the check-config endpoint and validates the response", async () => {
    let captured: unknown[] = [];
    const response = { status: "valid", config: { tree: "t.nwk", max_iter: 2 } };

    const bridge = createBridge(
      stubTransport({
        request: (...args) => {
          captured = args;

          return Promise.resolve(response);
        },
      }),
    );

    await expect(bridge.checkConfig({ command: "timetree", text: "tree: t.nwk\n" })).resolves.toStrictEqual(response);
    expect(captured).toStrictEqual(["check-config", { command: "timetree", text: "tree: t.nwk\n" }]);
  });

  test("checkConfig passes an invalid response through with its problems", async () => {
    const response = {
      status: "invalid",
      message: "invalid configuration: unknown field `x`",
      causes: [],
      problems: [{ code: "config::unknown-field", message: "unknown field `x`", span: { offset: 0, length: 1 } }],
      rendered: "x: 1",
    };

    const bridge = createBridge(stubTransport({ request: () => Promise.resolve(response) }));
    await expect(bridge.checkConfig({ command: "clock", text: "x: 1\n" })).resolves.toStrictEqual(response);
  });
});

describe("bridge commands and terminal events", () => {
  test("an ok terminal event resolves with the outcome and job events reach the callbacks", async () => {
    const events: JobEvent[] = [
      { type: "started", data: { job_id: "j1", command: "ancestral" } },
      { type: "progress", data: { stage: "read", fraction: 0.25, message: "reading" } },
      { type: "log", data: { level: "warn", message: "tip without date" } },
    ];

    let captured: unknown[] = [];

    const command: BridgeTransport["command"] = (name, config, options: TransportCommandOptions) => {
      captured = [name, config];

      for (const event of events) {
        options.onEvent(event);
      }

      return Promise.resolve({ status: "ok", job_id: "j1", result: OUTCOME });
    };

    const started: string[] = [];
    const progress: ProgressEvent[] = [];
    const logs: LogEvent[] = [];
    const bridge = createBridge(stubTransport({ command }));

    const result = await bridge.ancestral(
      { tree: "t.nwk", output_all: "out" },
      {
        onStarted: (jobId) => {
          started.push(jobId);
        },
        onProgress: (event) => {
          progress.push(event);
        },
        onLog: (event) => {
          logs.push(event);
        },
      },
    );

    expect(result).toStrictEqual(OUTCOME);
    expect(captured).toStrictEqual(["ancestral", { tree: "t.nwk", output_all: "out" }]);
    expect(started).toStrictEqual(["j1"]);
    expect(progress).toStrictEqual([{ stage: "read", fraction: 0.25, message: "reading" }]);
    expect(logs).toStrictEqual([{ level: "warn", message: "tip without date" }]);
  });

  test("an error terminal event rejects with a CommandError carrying the cause chain", async () => {
    const terminal = { status: "error", job_id: "j2", message: "When reading dates", causes: ["file not found"] };
    const bridge = createBridge(stubTransport({ command: () => Promise.resolve(terminal) }));
    const error: unknown = await bridge.clock({ tree: "t.nwk" }).catch((err: unknown) => err);

    expect(error).toBeInstanceOf(CommandError);
    expect(error).toMatchObject({ jobId: "j2", message: "When reading dates", causes: ["file not found"] });
  });

  test("a cancelled terminal event rejects with a CancelledError", async () => {
    const bridge = createBridge(
      stubTransport({ command: () => Promise.resolve({ status: "cancelled", job_id: "j3" }) }),
    );

    await expect(bridge.timetree({})).rejects.toBeInstanceOf(CancelledError);
  });

  test("a malformed terminal event rejects with a ZodError", async () => {
    const bridge = createBridge(stubTransport({ command: () => Promise.resolve({ status: "done" }) }));
    await expect(bridge.prune({ tree: "t.nwk" })).rejects.toBeInstanceOf(ZodError);
  });

  test("the abort signal is passed through to the transport", async () => {
    const controller = new AbortController();
    let captured: AbortSignal | undefined;

    const command: BridgeTransport["command"] = (_name, _config, options) => {
      captured = options.signal;

      return Promise.resolve({ status: "ok", job_id: "j4", result: { command: "optimize", output_files: [] } });
    };

    const bridge = createBridge(stubTransport({ command }));
    await bridge.optimize({ tree: "t.nwk" }, { signal: controller.signal });
    expect(captured).toBe(controller.signal);
  });

  test("cancelled error has a stable Error name", () => {
    const err = new CancelledError();
    expect(err).toBeInstanceOf(Error);
    expect(err.name).toBe("CancelledError");
  });
});

describe("bridge job event parser", () => {
  test("parseJobEvent accepts a terminal event", () => {
    const event = { type: "terminal", data: { status: "cancelled", job_id: "j" } };
    expect(parseJobEvent(event)).toStrictEqual(event);
  });

  test("parseJobEvent rejects an unknown event type", () => {
    expect(() => parseJobEvent({ type: "result", data: {} })).toThrow(ZodError);
  });

  test("parseJobEvent rejects a log level spelled in another case", () => {
    expect(() => parseJobEvent({ type: "log", data: { level: "Info", message: "hello" } })).toThrow(ZodError);
  });

  test("parseJobEvent rejects a progress event without a fraction", () => {
    expect(() => parseJobEvent({ type: "progress", data: { stage: "read", message: "half" } })).toThrow(ZodError);
  });
});
