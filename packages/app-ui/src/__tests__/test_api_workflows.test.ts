import { CancelledError, CommandError } from "@neherlab/app-contracts";
import { StreamEndedError } from "@neherlab/app-contracts/client";
import { describe, expect, test } from "vitest";

import { followRun, runCommand } from "../api/workflows";
import type { RunEvent } from "../results/progress";
import { eventStream, FakeServer, json, noDelay, RECORD, runEvent, type Route } from "./api_server";

const OUTCOME = { command: "clock", output_files: [{ path: "out/clock.nwk", kind: "nwk" }] };

const PROGRESS = [
  runEvent(0, "started", { job_id: "r1", command: "clock" }),
  runEvent(1, "progress", { stage: "clock", fraction: 0.5, message: "fitting" }),
];

function terminal(data: Record<string, unknown>) {
  return runEvent(2, "terminal", { job_id: "r1", ...data });
}

function server(events: Route, extra: Record<string, Route> = {}): FakeServer {
  return new FakeServer({
    "POST /api/runs": () => json(RECORD, 201),
    "GET /api/runs/r1/events?from=0": events,
    ...extra,
  });
}

describe("api workflows run a command", () => {
  test("a command creates a run, follows its events, and returns the outcome of the terminal event", async () => {
    const backend = server(() => eventStream([...PROGRESS, terminal({ status: "ok", result: OUTCOME })]));
    const started: string[] = [];
    const received: RunEvent[] = [];

    const outcome = await runCommand(
      backend.client(),
      { command: "clock", config: { tree: "t.nwk" }, title: "my clock" },
      {
        onStarted: (id) => {
          started.push(id);
        },
        onEvent: (event) => {
          received.push(event);
        },
        sleep: noDelay,
      },
    );

    expect(outcome).toStrictEqual(OUTCOME);
    expect(started).toStrictEqual(["r1"]);
    expect(received.map((event) => event.seq)).toStrictEqual([0, 1, 2]);
    expect(backend.sent).toStrictEqual([
      {
        key: "POST /api/runs",
        body: JSON.stringify({ command: "clock", config: { tree: "t.nwk" }, title: "my clock", defer_start: false }),
      },
      { key: "GET /api/runs/r1/events?from=0", body: "" },
    ]);
  });

  test.each([
    {
      status: "error",
      data: { status: "error", message: "When running clock", causes: ["no dates"] },
      type: CommandError,
      expected: { name: "CommandError", jobId: "r1", message: "When running clock", causes: ["no dates"] },
    },
    {
      status: "interrupted",
      data: { status: "interrupted" },
      type: CommandError,
      expected: {
        name: "CommandError",
        jobId: "r1",
        message: "the run was interrupted because the process that ran it stopped",
        causes: [],
      },
    },
    {
      status: "cancelled",
      data: { status: "cancelled" },
      type: CancelledError,
      expected: { name: "CancelledError", message: "Operation cancelled" },
    },
  ])("a $status terminal event rejects the command", async ({ data, type, expected }) => {
    const backend = server(() => eventStream([...PROGRESS, terminal(data)]));

    const error = await runCommand(backend.client(), { command: "clock", config: {} }, { sleep: noDelay }).catch(
      (failure: unknown) => failure,
    );

    expect(error).toBeInstanceOf(type);
    expect(error).toMatchObject(expected);
  });

  test("aborting the command asks the back end to cancel the run and rejects with the cancelled terminal event", async () => {
    const controller = new AbortController();
    let cancelled: () => void = () => undefined;

    const cancelRequested = new Promise<void>((resolve) => {
      cancelled = resolve;
    });

    const encoder = new TextEncoder();

    const backend = server(
      () =>
        new Response(
          new ReadableStream<Uint8Array>({
            async start(stream) {
              stream.enqueue(encoder.encode(PROGRESS.map((event) => `data: ${JSON.stringify(event)}\n\n`).join("")));
              await cancelRequested;
              stream.enqueue(encoder.encode(`data: ${JSON.stringify(terminal({ status: "cancelled" }))}\n\n`));
              stream.close();
            },
          }),
          { headers: { "Content-Type": "text/event-stream" } },
        ),
      {
        "POST /api/runs/r1/cancel": () => {
          cancelled();

          return json({ cancelled: true });
        },
      },
    );

    const command = runCommand(
      backend.client(),
      { command: "clock", config: {} },
      {
        signal: controller.signal,
        onEvent: (event) => {
          if (event.seq === 1) {
            controller.abort();
          }
        },
        sleep: noDelay,
      },
    );

    await expect(command).rejects.toBeInstanceOf(CancelledError);
    expect(backend.keys()).toStrictEqual([
      "POST /api/runs",
      "GET /api/runs/r1/events?from=0",
      "POST /api/runs/r1/cancel",
    ]);
  });

  test("a command aborted before it starts sends nothing", async () => {
    const backend = server(() => eventStream([]));

    await expect(
      runCommand(backend.client(), { command: "clock", config: {} }, { signal: AbortSignal.abort(), sleep: noDelay }),
    ).rejects.toBeInstanceOf(CancelledError);
    expect(backend.sent).toStrictEqual([]);
  });

  test("a run whose event stream stops delivering events without a terminal event rejects the command", async () => {
    const backend = server(() => eventStream(PROGRESS), { "GET /api/runs/r1/events?from=2": () => eventStream([]) });

    await expect(
      runCommand(backend.client(), { command: "clock", config: {} }, { sleep: noDelay }),
    ).rejects.toBeInstanceOf(StreamEndedError);
  });
});

describe("api workflows follow a run", () => {
  test("following a run from an offset returns its terminal event", async () => {
    const backend = new FakeServer({
      "GET /api/runs/r1/events?from=2": () => eventStream([terminal({ status: "ok", result: OUTCOME })]),
    });

    await expect(followRun(backend.client(), "r1", { from: 2, sleep: noDelay })).resolves.toStrictEqual({
      job_id: "r1",
      status: "ok",
      result: OUTCOME,
    });
  });

  test("following a run that is aborted before its terminal event rejects as cancelled", async () => {
    const controller = new AbortController();
    const backend = new FakeServer({ "GET /api/runs/r1/events?from=0": () => eventStream(PROGRESS) });

    await expect(
      followRun(backend.client(), "r1", {
        signal: controller.signal,
        onEvent: () => controller.abort(),
        sleep: noDelay,
      }),
    ).rejects.toBeInstanceOf(CancelledError);
  });
});
