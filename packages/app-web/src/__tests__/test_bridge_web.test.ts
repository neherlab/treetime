import { CancelledError } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";
import { ZodError } from "zod";

import { createWebBridge } from "../bridge-web";

interface StreamMessage {
  event: string;
  data: unknown;
}

const abortingFetch: typeof fetch = () => Promise.reject(new DOMException("Aborted", "AbortError"));

function eventStreamFetch(messages: StreamMessage[]): typeof fetch {
  const body = messages.map((message) => `event: ${message.event}\ndata: ${JSON.stringify(message.data)}\n\n`).join("");

  return () => Promise.resolve(new Response(body, { headers: { "Content-Type": "text/event-stream" } }));
}

function jsonFetch(body: unknown, init?: { status?: number; statusText?: string }): typeof fetch {
  return () => Promise.resolve(new Response(JSON.stringify(body), init));
}

describe("bridge_web query path", () => {
  test("version fetches and validates the result", async () => {
    const bridge = createWebBridge({ fetchFn: jsonFetch({ version: "9.9.9" }), apiBase: "/api" });
    await expect(bridge.version()).resolves.toStrictEqual({ version: "9.9.9" });
  });

  test("a non-ok response rejects with the status", async () => {
    const bridge = createWebBridge({ fetchFn: jsonFetch({}, { status: 500, statusText: "Server Error" }) });
    await expect(bridge.version()).rejects.toThrow("500");
  });
});

describe("bridge_web streaming command path", () => {
  test("progress and log events stream while the result resolves", async () => {
    const fetchFn = eventStreamFetch([
      { event: "progress", data: { stage: "read", fraction: 0.5, message: "reading" } },
      { event: "log", data: { level: "Info", message: "working" } },
      { event: "result", data: { model_name: "JC69" } },
    ]);

    const received: string[] = [];
    const bridge = createWebBridge({ fetchFn });

    const result = await bridge.ancestral(
      { tree: "t", outdir: "o" },
      {
        onProgress: (e) => {
          received.push(e.stage);
        },
      },
    );

    expect(result).toStrictEqual({ model_name: "JC69" });
    expect(received).toStrictEqual(["read"]);
  });

  test("a stream without a result event rejects", async () => {
    const fetchFn = eventStreamFetch([
      { event: "progress", data: { stage: "read", fraction: 0.5, message: "reading" } },
    ]);

    const bridge = createWebBridge({ fetchFn });
    await expect(bridge.ancestral({ tree: "t", outdir: "o" })).rejects.toThrow("no result received");
  });

  test("a non-ok command response rejects with the status", async () => {
    const bridge = createWebBridge({ fetchFn: jsonFetch({}, { status: 503, statusText: "Busy" }) });
    await expect(bridge.ancestral({ tree: "t", outdir: "o" })).rejects.toThrow("503");
  });

  test("a malformed result event rejects with a ZodError", async () => {
    const bridge = createWebBridge({ fetchFn: eventStreamFetch([{ event: "result", data: { wrong: true } }]) });
    await expect(bridge.ancestral({ tree: "t", outdir: "o" })).rejects.toBeInstanceOf(ZodError);
  });

  test("an abort surfaces as CancelledError", async () => {
    const bridge = createWebBridge({ fetchFn: abortingFetch });
    await expect(bridge.clock({ dates: "d", outdir: "o" })).rejects.toBeInstanceOf(CancelledError);
  });
});
