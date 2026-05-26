<<<<<<< HEAD
import type { TreeTimeBridge, CommandResult, VersionInfo } from "@neherlab/app-contracts";
=======
import { fetchEventSource } from "@microsoft/fetch-event-source";
import { CancelledError } from "@neherlab/app-contracts";
import type {
  CommandOptions,
  TreeTimeBridge,
  DatasetInfo,
  LogEvent,
  ProgressEvent,
  VersionInfo,
} from "@neherlab/app-contracts";
>>>>>>> ac719231 (feat(web): wire AbortController through bridge contract and UI)

const API_BASE = "/api";

async function getJson<T>(path: string): Promise<T> {
  const response = await fetch(`${API_BASE}/${path}`);
  if (!response.ok) {
    throw new Error(`GET ${path}: ${response.status} ${response.statusText}`);
  }
  return (await response.json()) as T;
}

<<<<<<< HEAD
async function postCommand(command: string, args: unknown): Promise<CommandResult> {
  const response = await fetch(`${API_BASE}/${command}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(args),
  });
<<<<<<< HEAD
=======
async function postSse<T>(command: string, args: unknown, options?: CommandOptions): Promise<T> {
  let result: T | undefined;

  try {
    await fetchEventSource(`${API_BASE}/${command}`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(args),
      signal: options?.signal,
      onmessage(msg) {
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
>>>>>>> ac719231 (feat(web): wire AbortController through bridge contract and UI)

  if (!response.ok) {
    const text = await response.text();
    return { status: "error", error: text };
  }

  return response.json() as Promise<CommandResult>;
=======
  if (!response.ok) {
    return { status: "error", error: `${response.status} ${response.statusText}` };
  }
  return (await response.json()) as CommandResult;
>>>>>>> 33bee034 (feat: add end-to-end version info across all layers)
}

export function createWebBridge(): TreeTimeBridge {
  return {
    version: () => getJson<VersionInfo>("version"),
<<<<<<< HEAD
    ancestral: (args) => postCommand("ancestral", args),
    clock: (args) => postCommand("clock", args),
    timetree: (args) => postCommand("timetree", args),
    mugration: (args) => postCommand("mugration", args),
    optimize: (args) => postCommand("optimize", args),
    prune: (args) => postCommand("prune", args),
=======
    datasets: () => getJson<DatasetInfo[]>("datasets"),
    ancestral: (args, options) => postSse("ancestral", args, options),
    clock: (args, options) => postSse("clock", args, options),
    timetree: (args, options) => postSse("timetree", args, options),
    mugration: (args, options) => postSse("mugration", args, options),
    optimize: (args, options) => postSse("optimize", args, options),
    prune: (args, options) => postSse("prune", args, options),
>>>>>>> ac719231 (feat(web): wire AbortController through bridge contract and UI)
  };
}
