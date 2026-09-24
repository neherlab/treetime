import { CancelledError, CommandError } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";
import { ZodError } from "zod";

import { createWebBridge } from "../bridge-web";

interface StreamMessage {
  event: string;
  data: unknown;
}

interface Call {
  url: string;
  body: BodyInit | null | undefined;
}

const OUTCOME = { command: "ancestral", output_files: ["out/ancestral.nwk"] };

function eventStreamFetch(messages: StreamMessage[], calls: Call[] = []): typeof fetch {
  return (input, init) => {
    calls.push({ url: urlOf(input), body: init?.body });

    return Promise.resolve(new Response(sse(messages), { headers: { "Content-Type": "text/event-stream" } }));
  };
}

function sse(messages: StreamMessage[]): string {
  return messages.map((message) => `event: ${message.event}\ndata: ${JSON.stringify(message.data)}\n\n`).join("");
}

function urlOf(input: RequestInfo | URL): string {
  if (input instanceof Request) {
    return input.url;
  }

  return input instanceof URL ? input.href : input;
}

function jsonFetch(body: unknown, init?: { status?: number; statusText?: string }): typeof fetch {
  return () => Promise.resolve(new Response(JSON.stringify(body), init));
}

describe("bridge_web query and request paths", () => {
  test("version fetches and validates the result", async () => {
    const bridge = createWebBridge({ fetchFn: jsonFetch({ version: "9.9.9" }), apiBase: "/api" });
    await expect(bridge.version()).resolves.toStrictEqual({ version: "9.9.9" });
  });

  test("a non-ok response rejects with the status", async () => {
    const bridge = createWebBridge({ fetchFn: jsonFetch({}, { status: 500, statusText: "Server Error" }) });
    await expect(bridge.version()).rejects.toThrow("500");
  });

  test("checkConfig posts the request body", async () => {
    const calls: Call[] = [];

    const fetchFn: typeof fetch = (input, init) => {
      calls.push({ url: urlOf(input), body: init?.body });

      return Promise.resolve(new Response(JSON.stringify({ status: "valid", config: { tree: "t.nwk" } })));
    };

    const bridge = createWebBridge({ fetchFn });
    await expect(bridge.checkConfig({ command: "prune", text: "tree: t.nwk" })).resolves.toStrictEqual({
      status: "valid",
      config: { tree: "t.nwk" },
    });
    expect(calls).toStrictEqual([
      { url: "/api/check-config", body: JSON.stringify({ command: "prune", text: "tree: t.nwk" }) },
    ]);
  });
});

describe("bridge_web streaming command path", () => {
  test("job events stream while the terminal event resolves the outcome", async () => {
    const calls: Call[] = [];

    const fetchFn = eventStreamFetch(
      [
        { event: "started", data: { job_id: "j1", command: "ancestral" } },
        { event: "progress", data: { stage: "read", fraction: 0.5, message: "reading" } },
        { event: "log", data: { level: "info", message: "working" } },
        { event: "terminal", data: { status: "ok", job_id: "j1", result: OUTCOME } },
      ],
      calls,
    );

    const received: string[] = [];
    const bridge = createWebBridge({ fetchFn });

    const result = await bridge.ancestral(
      { tree: "t.nwk" },
      {
        onStarted: (jobId) => {
          received.push(`started ${jobId}`);
        },
        onProgress: (event) => {
          received.push(event.stage);
        },
      },
    );

    expect(result).toStrictEqual(OUTCOME);
    expect(received).toStrictEqual(["started j1", "read"]);
    expect(calls).toStrictEqual([{ url: "/api/ancestral", body: JSON.stringify({ tree: "t.nwk" }) }]);
  });

  test("an error terminal event rejects with a CommandError", async () => {
    const fetchFn = eventStreamFetch([
      { event: "started", data: { job_id: "j2", command: "clock" } },
      { event: "terminal", data: { status: "error", job_id: "j2", message: "bad", causes: ["cause"] } },
    ]);

    const bridge = createWebBridge({ fetchFn });
    await expect(bridge.clock({ tree: "t.nwk" })).rejects.toBeInstanceOf(CommandError);
  });

  test("a stream without a terminal event rejects", async () => {
    const fetchFn = eventStreamFetch([{ event: "started", data: { job_id: "j3", command: "ancestral" } }]);
    const bridge = createWebBridge({ fetchFn });
    await expect(bridge.ancestral({ tree: "t.nwk" })).rejects.toThrow("ended without a terminal event");
  });

  test("a non-ok command response rejects with the status", async () => {
    const bridge = createWebBridge({ fetchFn: jsonFetch({}, { status: 503, statusText: "Busy" }) });
    await expect(bridge.ancestral({ tree: "t.nwk" })).rejects.toThrow("503");
  });

  test("a malformed event rejects with a ZodError", async () => {
    const bridge = createWebBridge({ fetchFn: eventStreamFetch([{ event: "result", data: { wrong: true } }]) });
    await expect(bridge.ancestral({ tree: "t.nwk" })).rejects.toBeInstanceOf(ZodError);
  });

  test("aborting requests cancellation of the job and resolves with its cancelled terminal event", async () => {
    const controller = new AbortController();
    const calls: string[] = [];
    let release: () => void = () => undefined;

    const released = new Promise<void>((resolve) => {
      release = resolve;
    });

    const fetchFn: typeof fetch = (input) => {
      const url = urlOf(input);
      calls.push(url);

      if (url.endsWith("/cancel")) {
        release();

        return Promise.resolve(new Response(JSON.stringify({ cancelled: true })));
      }

      const encoder = new TextEncoder();

      const body = new ReadableStream<Uint8Array>({
        async start(stream) {
          stream.enqueue(encoder.encode(sse([{ event: "started", data: { job_id: "j4", command: "timetree" } }])));
          controller.abort();
          await released;
          stream.enqueue(encoder.encode(sse([{ event: "terminal", data: { status: "cancelled", job_id: "j4" } }])));
          stream.close();
        },
      });

      return Promise.resolve(new Response(body, { headers: { "Content-Type": "text/event-stream" } }));
    };

    const bridge = createWebBridge({ fetchFn });
    await expect(bridge.timetree({}, { signal: controller.signal })).rejects.toBeInstanceOf(CancelledError);
    expect(calls).toStrictEqual(["/api/timetree", "/api/jobs/j4/cancel"]);
  });
});
