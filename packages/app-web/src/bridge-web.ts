import { fetchEventSource } from "@microsoft/fetch-event-source";
import {
  CancelledError,
  createBridge,
  parseLogEvent,
  parseProgressEvent,
  type BridgeTransport,
  type CommandOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";

export interface WebBridgeDeps {
  fetchFn?: typeof fetch;
  fetchEventSourceFn?: typeof fetchEventSource;
  apiBase?: string;
}

export function createWebBridge(deps: WebBridgeDeps = {}): TreeTimeBridge {
  const fetchFn = deps.fetchFn ?? globalThis.fetch.bind(globalThis);
  const fetchEventSourceFn = deps.fetchEventSourceFn ?? fetchEventSource;
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
      await fetchEventSourceFn(`${apiBase}/${command}`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(args),
        signal: options?.signal ?? null,
        onmessage(msg) {
          if (debug) console.debug("[TreeTime]", JSON.stringify(msg));
          if (msg.event === "progress") {
            options?.onProgress?.(parseProgressEvent(JSON.parse(msg.data)));
          } else if (msg.event === "log") {
            logToConsole(parseLogEvent(JSON.parse(msg.data)));
          } else if (msg.event === "result") {
            result = JSON.parse(msg.data);
            received = true;
          }
        },
        onerror(err) {
          throw err;
        },
        openWhenHidden: true,
      });
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
  return env.TREETIME_DEBUG_FETCH === "true" || (env.DEV && env.TREETIME_DEBUG_FETCH !== "false");
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
