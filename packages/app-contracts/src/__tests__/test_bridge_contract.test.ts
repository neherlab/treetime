import { describe, expect, test } from "vitest";
import { ZodError } from "zod";

import {
  CancelledError,
  createBridge,
  parseLogEvent,
  parseProgressEvent,
  type BridgeTransport,
  type CommandOptions,
  type ProgressEvent,
} from "../index";

function stubTransport(overrides: Partial<BridgeTransport>): BridgeTransport {
  return {
    query: overrides.query ?? (() => Promise.reject(new Error("query not stubbed"))),
    command: overrides.command ?? (() => Promise.reject(new Error("command not stubbed"))),
  };
}

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

  test("a malformed command result rejects with a ZodError", async () => {
    const bridge = createBridge(stubTransport({ command: () => Promise.resolve({ wrong: true }) }));
    await expect(bridge.ancestral({ tree: "t", outdir: "o" })).rejects.toBeInstanceOf(ZodError);
  });

  test("a malformed query result rejects with a ZodError", async () => {
    const bridge = createBridge(stubTransport({ query: () => Promise.resolve({ version: 1 }) }));
    await expect(bridge.version()).rejects.toBeInstanceOf(ZodError);
  });
});

describe("bridge streaming and cancellation", () => {
  test("command forwards progress events and returns the validated result", async () => {
    const progress: ProgressEvent[] = [
      { stage: "read", fraction: 0.25, message: "reading" },
      { stage: "infer", fraction: 1, message: "done" },
    ];

    const command: BridgeTransport["command"] = (_endpoint, _args, options?: CommandOptions) => {
      for (const event of progress) {
        options?.onProgress?.(event);
      }
      return Promise.resolve({ model_name: "JC69" });
    };

    const received: ProgressEvent[] = [];
    const bridge = createBridge(stubTransport({ command }));
    const result = await bridge.ancestral({ tree: "t", outdir: "o" }, { onProgress: (e) => received.push(e) });

    expect(result).toStrictEqual({ model_name: "JC69" });
    expect(received).toStrictEqual(progress);
  });

  test("a CancelledError from the transport propagates unchanged", async () => {
    const bridge = createBridge(stubTransport({ command: () => Promise.reject(new CancelledError()) }));
    await expect(bridge.clock({ dates: "d", outdir: "o" })).rejects.toBeInstanceOf(CancelledError);
  });

  test("CancelledError is an Error with a stable name", () => {
    const err = new CancelledError();
    expect(err).toBeInstanceOf(Error);
    expect(err.name).toBe("CancelledError");
  });

  test("the abort signal is passed through to the transport", async () => {
    const controller = new AbortController();
    let captured: Parameters<BridgeTransport["command"]> | undefined;
    const command: BridgeTransport["command"] = (...args) => {
      captured = args;
      return Promise.resolve({});
    };
    const bridge = createBridge(stubTransport({ command }));

    await bridge.optimize({ tree: "t", outdir: "o" }, { signal: controller.signal });

    expect(captured).toStrictEqual(["optimize", { tree: "t", outdir: "o" }, { signal: controller.signal }]);
  });
});

describe("bridge streaming event parsers", () => {
  test("parseProgressEvent accepts a well-formed event", () => {
    const event = { stage: "read", fraction: 0.5, message: "half" };
    expect(parseProgressEvent(event)).toStrictEqual(event);
  });

  test("parseProgressEvent rejects a missing fraction", () => {
    expect(() => parseProgressEvent({ stage: "read", message: "half" })).toThrow(ZodError);
  });

  test("parseLogEvent accepts a well-formed event", () => {
    const event = { level: "Info", message: "hello" };
    expect(parseLogEvent(event)).toStrictEqual(event);
  });

  test("parseLogEvent rejects an unknown level", () => {
    expect(() => parseLogEvent({ level: "Fatal", message: "hello" })).toThrow(ZodError);
  });
});
