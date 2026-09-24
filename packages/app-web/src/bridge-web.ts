import {
  CancelledError,
  createBridge,
  parseLogEvent,
  parseProgressEvent,
  type BridgeTransport,
  type CommandOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";
import { EventSourceParserStream } from "eventsource-parser/stream";

export interface WebBridgeDeps {
  fetchFn?: typeof fetch;
  apiBase?: string;
}

export function createWebBridge(deps: WebBridgeDeps = {}): TreeTimeBridge {
  const fetchFn = deps.fetchFn ?? globalThis.fetch.bind(globalThis);
  const apiBase = deps.apiBase ?? "/api";

  const debug = readDebugFlag();

  async function getJson(path: string): Promise<unknown> {
    if (debug) console.debug("[TreeTime] GET", path);
    const response = await fetchFn(`${apiBase}/${path}`);

    if (!response.ok) {
      throw new Error(`GET ${path}: ${response.status} ${response.statusText}`);
    }

    const data: unknown = await response.json();

    if (debug) console.debug("[TreeTime] GET", path, JSON.stringify(data));

    return data;
  }

  async function postSse(command: string, args: unknown, options?: CommandOptions): Promise<unknown> {
    let result: unknown;
    let received = false;

    try {
      const response = await fetchFn(`${apiBase}/${command}`, {
        method: "POST",
        headers: { "Content-Type": "application/json", Accept: "text/event-stream" },
        body: JSON.stringify(args),
        signal: options?.signal ?? null,
      });

      if (!response.ok || response.body === null) {
        throw new Error(`POST ${command}: ${response.status} ${response.statusText}`);
      }

      const events = response.body.pipeThrough(new TextDecoderStream()).pipeThrough(new EventSourceParserStream());

      for await (const message of events) {
        if (debug) console.debug("[TreeTime]", JSON.stringify(message));

        if (message.event === "progress") {
          options?.onProgress?.(parseProgressEvent(JSON.parse(message.data)));
        } else if (message.event === "log") {
          logToConsole(parseLogEvent(JSON.parse(message.data)));
        } else if (message.event === "result") {
          result = JSON.parse(message.data);
          received = true;
        }
      }
    } catch (err: unknown) {
      if (err instanceof DOMException && err.name === "AbortError") {
        throw new CancelledError();
      }

      throw err;
    }

    if (!received) {
      throw new Error(`${command}: no result received`);
    }

    return result;
  }

  const transport: BridgeTransport = {
    query: (endpoint) => getJson(endpoint),
    command: (endpoint, args, options) => postSse(endpoint, args, options),
  };

  return createBridge(transport);
}

function readDebugFlag(): boolean {
  const env = import.meta.env;

  return env.VITE_TREETIME_DEBUG_FETCH === "true" || (env.DEV && env.VITE_TREETIME_DEBUG_FETCH !== "false");
}

function logToConsole(log: { level: string; message: string }): void {
  switch (log.level) {
    case "Error":
      console.error(`[TreeTime] ${log.message}`);
      break;
    case "Warn":
      console.warn(`[TreeTime] ${log.message}`);
      break;
    default:
      console.log(`[TreeTime] [${log.level}] ${log.message}`);
      break;
  }
}
