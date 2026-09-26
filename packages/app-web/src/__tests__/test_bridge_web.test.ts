import {
  BridgeError,
  CancelledError,
  CommandError,
  type TreeTimeBridge,
  zOperationRequest,
} from "@neherlab/app-contracts";
import { ApiError, createApiClient } from "@neherlab/app-contracts/client";
import { describe, expect, test } from "vitest";
import { ZodError, z } from "zod";

import { createWebBridge } from "../bridge-web";

interface Call {
  method: string;
  url: string;
  body: BodyInit | null | undefined;
}

type Route = (call: Call) => Response;

interface TestBridgeDeps {
  fetchFn: typeof fetch;
  saveBlob?: (blob: Blob, name: string) => void;
}

const BASE_URL = "http://treetime.test";

const RECORD = {
  id: "r1",
  title: "clock",
  command: "clock",
  config: { tree: "t.nwk" },
  status: "running",
  pinned: false,
  created_at: "2026-09-25T10:00:00Z",
  started_at: null,
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

const OUTCOME = { command: "clock", output_files: [{ path: "/runs/r1/out/clock.nwk", kind: "nwk" }] };

function webBridge(deps: TestBridgeDeps): TreeTimeBridge {
  return createWebBridge({ ...deps, client: createApiClient({ baseUrl: BASE_URL, fetch: deps.fetchFn }) });
}

function routes(table: Record<string, Route>, calls: Call[] = []): typeof fetch {
  return async (input, init) => {
    const call =
      input instanceof Request
        ? { method: input.method, url: urlOf(input), body: await input.blob() }
        : { method: init?.method ?? "GET", url: urlOf(input), body: init?.body };

    calls.push(call);
    const route = table[routeKey(call)];

    if (route === undefined) {
      return new Response(JSON.stringify({ code: "x", message: "no route" }), { status: 404 });
    }

    if (init?.signal?.aborted === true) {
      throw new DOMException("aborted", "AbortError");
    }

    return route(call);
  };
}

function json(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), { status, headers: { "Content-Type": "application/json" } });
}

function sse(events: Array<{ type: string; data: unknown }>): Response {
  const text = events
    .map(
      (event, seq) =>
        `event: ${event.type}\nid: ${seq}\ndata: ${JSON.stringify({ seq, time: "t", type: event.type, data: event.data })}\n\n`,
    )
    .join("");

  return new Response(text, { headers: { "Content-Type": "text/event-stream" } });
}

function routeKey(call: Call): string {
  if (call.url === "/api/operations") {
    return `call ${zOperationRequest.parse(JSON.parse(z.string().parse(call.body))).operation}`;
  }

  return `${call.method} ${call.url}`;
}

function urlOf(input: RequestInfo | URL): string {
  if (input instanceof Request) {
    const url = new URL(input.url);

    return `${url.pathname}${url.search}`;
  }

  return input instanceof URL ? input.href : input;
}

describe("bridge_web queries and requests", () => {
  test("version fetches and validates the result", async () => {
    const bridge = webBridge({ fetchFn: routes({ "call version": () => json({ version: "9.9.9" }) }) });
    await expect(bridge.version()).resolves.toStrictEqual({ version: "9.9.9" });
  });

  test("a failed request rejects with the server's typed error, its causes and the status", async () => {
    const response = { code: "not_found", message: "When reading run `r9`", causes: ["no run with id `r9`"] };
    const bridge = webBridge({ fetchFn: routes({ "call get-run": () => json(response, 404) }) });

    const error = await bridge.getRun("r9").catch((failure: unknown) => failure);

    expect(error).toBeInstanceOf(BridgeError);
    expect(error).toMatchObject({
      message: "get-run: 404: When reading run `r9`: no run with id `r9`",
      response,
    });
  });

  test("a failed request without a typed error body rejects as an internal error with the body text", async () => {
    const bridge = webBridge({
      fetchFn: routes({ "call get-run": () => new Response("Bad Gateway", { status: 502 }) }),
    });

    const error = await bridge.getRun("r9").catch((failure: unknown) => failure);

    expect(error).toMatchObject({
      message: "get-run: 502: Bad Gateway",
      response: { code: "internal_error", message: "Bad Gateway", causes: [] },
    });
  });

  test("a malformed result rejects with a ZodError", async () => {
    const bridge = webBridge({ fetchFn: routes({ "call version": () => json({ version: 1 }) }) });
    await expect(bridge.version()).rejects.toBeInstanceOf(ZodError);
  });

  test("checkConfig posts the check-config operation request", async () => {
    const calls: Call[] = [];

    const fetchFn = routes(
      {
        "call check-config": () =>
          json({
            status: "valid",
            command: "prune",
            config: { tree: "t.nwk" },
            code: {
              command_line: [],
              command_line_text: "",
              yaml: [],
              yaml_text: "",
              config_file: "run.yaml",
              config_command: "treetime timetree --config run.yaml",
            },
            checks: [],
          }),
      },
      calls,
    );

    const bridge = webBridge({ fetchFn });
    await bridge.checkConfig({ command: "prune", text: "tree: t.nwk" });
    expect(calls).toStrictEqual([
      {
        method: "POST",
        url: "/api/operations",
        body: JSON.stringify({
          operation: "check-config",
          args: { request: { command: "prune", text: "tree: t.nwk" } },
        }),
      },
    ]);
  });

  test("runConfig posts the run-config operation request", async () => {
    const calls: Call[] = [];

    const fetchFn = routes(
      {
        "call run-config": () =>
          json({
            status: "valid",
            config: { tree: "t.nwk" },
            code: {
              command_line: [],
              command_line_text: "",
              yaml: [],
              yaml_text: "",
              config_file: "run.yaml",
              config_command: "treetime timetree --config run.yaml",
            },
          }),
      },
      calls,
    );

    const bridge = webBridge({ fetchFn });
    await bridge.runConfig({ command: "prune", config: { tree: "t.nwk" } });
    expect(calls).toStrictEqual([
      {
        method: "POST",
        url: "/api/operations",
        body: JSON.stringify({
          operation: "run-config",
          args: { request: { command: "prune", config: { tree: "t.nwk" } } },
        }),
      },
    ]);
  });

  test("run operations send their operation requests", async () => {
    const calls: Call[] = [];
    const summary = { ...RECORD, status: "ok" };

    const fetchFn = routes(
      {
        "call update-run": () => json(summary),
        "call cancel-run": () => json({ cancelled: false }),
        "call delete-run": () => json(null),
        "call restore-run": () => json(summary),
        "call start-run": () => json(RECORD),
      },
      calls,
    );

    const bridge = webBridge({ fetchFn });
    await bridge.updateRun("r1", { pinned: true });
    await expect(bridge.cancelRun("r1")).resolves.toBe(false);
    await bridge.deleteRun("r1");
    await bridge.restoreRun("r1");
    await bridge.startRun("r1", { config: { tree: "/runs/r1/inputs/t.nwk" } });
    expect(calls.map((call) => call.body)).toStrictEqual([
      JSON.stringify({ operation: "update-run", args: { id: "r1", request: { pinned: true } } }),
      JSON.stringify({ operation: "cancel-run", args: { id: "r1" } }),
      JSON.stringify({ operation: "delete-run", args: { id: "r1" } }),
      JSON.stringify({ operation: "restore-run", args: { id: "r1" } }),
      JSON.stringify({
        operation: "start-run",
        args: { id: "r1", request: { config: { tree: "/runs/r1/inputs/t.nwk" } } },
      }),
    ]);
  });

  test("uploadInput puts the file bytes", async () => {
    const calls: Call[] = [];
    const uploaded = { name: "a b.nwk", path: "/runs/r1/inputs/a b.nwk", size: 6, sha256: "x" };
    const fetchFn = routes({ "PUT /api/runs/r1/inputs/a%20b.nwk": () => json(uploaded) }, calls);
    const bridge = webBridge({ fetchFn });
    await expect(bridge.uploadInput("r1", "a b.nwk", new Blob(["(A,B);"]))).resolves.toStrictEqual(uploaded);
    const body = calls[0]?.body;
    expect(body).toBeInstanceOf(Blob);
    expect(body instanceof Blob ? await body.text() : undefined).toBe("(A,B);");
  });

  test("a rejected upload rejects with the server's ApiError", async () => {
    const response = { code: "upload_too_large", message: "The file is too large", causes: [] };
    const fetchFn = routes({ "PUT /api/runs/r1/inputs/t.nwk": () => json(response, 413) });
    const bridge = webBridge({ fetchFn });
    const error = await bridge.uploadInput("r1", "t.nwk", new Blob(["(A,B);"])).catch((failure: unknown) => failure);
    expect(error).toBeInstanceOf(ApiError);
    expect(error).toMatchObject({ status: 413, response });
  });

  test("saveRunFile hands the file to the browser download under the given name", async () => {
    const saved: Array<{ name: string; blob: Blob }> = [];
    const fetchFn = routes({ "GET /api/runs/r1/file?path=out%2Fclock.nwk": () => new Response("(A,B);") });

    const bridge = webBridge({
      fetchFn,
      saveBlob: (blob, name) => {
        saved.push({ name, blob });
      },
    });

    await expect(bridge.saveRunFile("r1", "out/clock.nwk", "clock.nwk")).resolves.toBe(true);
    const contents = await Promise.all(saved.map(async ({ name, blob }) => [name, await blob.text()]));
    expect(contents).toStrictEqual([["clock.nwk", "(A,B);"]]);
  });

  test("saveRunArchive downloads the archive of the run", async () => {
    const names: string[] = [];
    const fetchFn = routes({ "GET /api/runs/r1/archive": () => new Response("PK") });

    const bridge = webBridge({
      fetchFn,
      saveBlob: (_blob, name) => {
        names.push(name);
      },
    });

    await expect(bridge.saveRunArchive("r1", "clock.zip")).resolves.toBe(true);
    expect(names).toStrictEqual(["clock.zip"]);
  });

  test("a failed download rejects without saving", async () => {
    const names: string[] = [];

    const fetchFn = routes({
      "GET /api/runs/r9/archive": () => json({ code: "not_found", message: "no run with id `r9`", causes: [] }, 404),
    });

    const bridge = webBridge({
      fetchFn,
      saveBlob: (_blob, name) => {
        names.push(name);
      },
    });

    await expect(bridge.saveRunArchive("r9", "r9.zip")).rejects.toThrow(
      "GET runs/r9/archive: 404: no run with id `r9`",
    );
    expect(names).toStrictEqual([]);
  });
});

describe("bridge_web run events", () => {
  test("followRun reads the event stream from the given offset", async () => {
    const calls: Call[] = [];

    const fetchFn = routes(
      {
        "GET /api/runs/r1/events?from=0": () =>
          sse([
            { type: "started", data: { job_id: "r1", command: "clock" } },
            { type: "log", data: { level: "info", message: "hello" } },
            { type: "terminal", data: { status: "cancelled", job_id: "r1" } },
          ]),
      },
      calls,
    );

    const bridge = webBridge({ fetchFn });
    const seen: string[] = [];

    const terminal = await bridge.followRun("r1", {
      onEvent: (event) => {
        seen.push(event.type);
      },
    });

    expect(terminal).toStrictEqual({ status: "cancelled", job_id: "r1" });
    expect(seen).toStrictEqual(["started", "log", "terminal"]);
  });

  test("a command creates a run and resolves with the outcome of its terminal event", async () => {
    const fetchFn = routes({
      "call create-run": () => json(RECORD),
      "GET /api/runs/r1/events?from=0": () =>
        sse([
          { type: "started", data: { job_id: "r1", command: "clock" } },
          { type: "progress", data: { stage: "read", fraction: 0.5, message: "" } },
          { type: "terminal", data: { status: "ok", job_id: "r1", result: OUTCOME } },
        ]),
    });

    const bridge = webBridge({ fetchFn });
    const fractions: number[] = [];
    await expect(
      bridge.clock(
        { tree: "t.nwk" },
        {
          onProgress: (event) => {
            fractions.push(event.fraction);
          },
        },
      ),
    ).resolves.toStrictEqual(OUTCOME);
    expect(fractions).toStrictEqual([0.5]);
  });

  test("an error terminal event rejects with a CommandError", async () => {
    const fetchFn = routes({
      "call create-run": () => json(RECORD),
      "GET /api/runs/r1/events?from=0": () =>
        sse([{ type: "terminal", data: { status: "error", job_id: "r1", message: "bad", causes: [] } }]),
    });

    const bridge = webBridge({ fetchFn });
    await expect(bridge.clock({ tree: "t.nwk" })).rejects.toBeInstanceOf(CommandError);
  });

  test("an aborted stream rejects with a CancelledError", async () => {
    const controller = new AbortController();
    controller.abort();
    const fetchFn = routes({ "GET /api/runs/r1/events?from=2": () => sse([]) });
    const bridge = webBridge({ fetchFn });
    await expect(bridge.followRun("r1", { from: 2, signal: controller.signal })).rejects.toBeInstanceOf(CancelledError);
  });
});
