import type { fetchEventSource } from "@microsoft/fetch-event-source";
import { CancelledError } from "@neherlab/app-contracts";
import { describe, expect, test } from "vitest";
import { ZodError } from "zod";

import { createWebBridge } from "../bridge-web";

type FetchEventSource = typeof fetchEventSource;

type FetchEventSourceInit = Parameters<FetchEventSource>[1];

interface StreamMessage {
  event: string;
  data: unknown;
}

const rejectingEventSource: FetchEventSource = () => Promise.reject(new DOMException("Aborted", "AbortError"));

function scriptedEventSource(messages: StreamMessage[]): FetchEventSource {
  return (_input, init: FetchEventSourceInit) => {
    for (const message of messages) {
      init.onmessage?.({ id: "", event: message.event, data: JSON.stringify(message.data) });
    }

    return Promise.resolve();
  };
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
    const eventSource = scriptedEventSource([
      { event: "progress", data: { stage: "read", fraction: 0.5, message: "reading" } },
      { event: "log", data: { level: "Info", message: "working" } },
      { event: "result", data: { model_name: "JC69" } },
    ]);

    const received: string[] = [];
    const bridge = createWebBridge({ fetchEventSourceFn: eventSource });
    const result = await bridge.ancestral({ tree: "t", outdir: "o" }, { onProgress: (e) => received.push(e.stage) });

    expect(result).toStrictEqual({ model_name: "JC69" });
    expect(received).toStrictEqual(["read"]);
  });

  test("a stream without a result event rejects", async () => {
    const eventSource = scriptedEventSource([
      { event: "progress", data: { stage: "read", fraction: 0.5, message: "reading" } },
    ]);

    const bridge = createWebBridge({ fetchEventSourceFn: eventSource });
    await expect(bridge.ancestral({ tree: "t", outdir: "o" })).rejects.toThrow("no result received");
  });

  test("a malformed result event rejects with a ZodError", async () => {
    const eventSource = scriptedEventSource([{ event: "result", data: { wrong: true } }]);
    const bridge = createWebBridge({ fetchEventSourceFn: eventSource });
    await expect(bridge.ancestral({ tree: "t", outdir: "o" })).rejects.toBeInstanceOf(ZodError);
  });

  test("an abort surfaces as CancelledError", async () => {
    const bridge = createWebBridge({ fetchEventSourceFn: rejectingEventSource });
    await expect(bridge.clock({ dates: "d", outdir: "o" })).rejects.toBeInstanceOf(CancelledError);
  });
});
