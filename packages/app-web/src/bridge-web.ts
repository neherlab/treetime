import { fetchEventSource } from "@microsoft/fetch-event-source";
import {
  CancelledError,
  createBridge,
  type BridgeTransport,
  type CommandOptions,
  type LogEvent,
  type ProgressEvent,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";

const API_BASE = "/api";

const DEBUG_FETCH =
  import.meta.env.TREETIME_DEBUG_FETCH === "true" ||
  (import.meta.env.DEV && import.meta.env.TREETIME_DEBUG_FETCH !== "false");

async function getJson<T>(path: string): Promise<T> {
  if (DEBUG_FETCH) console.debug("[TreeTime] GET", path);
  const response = await fetch(`${API_BASE}/${path}`);
  if (!response.ok) {
    throw new Error(`GET ${path}: ${response.status} ${response.statusText}`);
  }
  const data = (await response.json()) as T;
  if (DEBUG_FETCH) console.debug("[TreeTime] GET", path, JSON.stringify(data));
  return data;
}

async function postSse<T>(command: string, args: unknown, options?: CommandOptions): Promise<T> {
  let result: T | undefined;

  try {
    await fetchEventSource(`${API_BASE}/${command}`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(args),
      signal: options?.signal,
      onmessage(msg) {
        if (DEBUG_FETCH) console.debug("[TreeTime]", JSON.stringify(msg));
        if (msg.event === "progress") {
          options?.onProgress?.(JSON.parse(msg.data) as ProgressEvent);
        } else if (msg.event === "log") {
          const log = JSON.parse(msg.data) as LogEvent;
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
        } else if (msg.event === "result") {
          result = JSON.parse(msg.data) as T;
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

  if (result === undefined) {
    throw new Error(`${command}: no result received`);
  }

  return result;
}

function createWebTransport(): BridgeTransport {
  return {
    query: <T>(endpoint: string) => getJson<T>(endpoint),
    command: <T>(endpoint: string, args: unknown, options?: CommandOptions) => postSse<T>(endpoint, args, options),
  };
}

export function createWebBridge(): TreeTimeBridge {
  return createBridge(createWebTransport());
}
