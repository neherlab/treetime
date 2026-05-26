<<<<<<< HEAD
<<<<<<< HEAD
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
=======
import { fetchEventSource } from "@microsoft/fetch-event-source";
<<<<<<< HEAD
import type { TreeTimeBridge, ProgressEvent, VersionInfo } from "@neherlab/app-contracts";
>>>>>>> 3e5b08aa (feat(app): wire real web bridge with SSE progress transport)
=======
import type { ErrorResponse, TreeTimeBridge, VersionInfo } from "@neherlab/app-contracts";
>>>>>>> c2b9da5e (feat: add per-command result types and hooks across TypeScript layer)
=======
import type { TreeTimeBridge, DatasetInfo, ProgressEvent, VersionInfo } from "@neherlab/app-contracts";
>>>>>>> b8625b9a (feat(app): wire datasets through bridge contract and implementations)

const API_BASE = "/api";

async function getJson<T>(path: string): Promise<T> {
  const response = await fetch(`${API_BASE}/${path}`);
  if (!response.ok) {
    throw new Error(`GET ${path}: ${response.status} ${response.statusText}`);
  }
  return (await response.json()) as T;
}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
async function postCommand(command: string, args: unknown): Promise<CommandResult> {
=======
async function postCommand<T>(command: string, args: unknown): Promise<T> {
>>>>>>> c2b9da5e (feat: add per-command result types and hooks across TypeScript layer)
  const response = await fetch(`${API_BASE}/${command}`, {
=======
async function postSse<T>(command: string, args: unknown, onProgress?: (event: ProgressEvent) => void): Promise<T> {
  let result: T | undefined;

  await fetchEventSource(`${API_BASE}/${command}`, {
>>>>>>> 3e5b08aa (feat(app): wire real web bridge with SSE progress transport)
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(args),
    onmessage(msg) {
      if (msg.event === "progress") {
        onProgress?.(JSON.parse(msg.data) as ProgressEvent);
      } else if (msg.event === "result") {
        result = JSON.parse(msg.data) as T;
      }
    },
    onerror(err) {
      throw err;
    },
    openWhenHidden: true,
  });
<<<<<<< HEAD
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
    const body = (await response.json()) as ErrorResponse;
    throw new Error(`${body.code}: ${body.message}`);
  }
<<<<<<< HEAD
  return (await response.json()) as CommandResult;
>>>>>>> 33bee034 (feat: add end-to-end version info across all layers)
=======

  if (result === undefined) {
    throw new Error(`${command}: no result received`);
  }

  return result;
>>>>>>> 3e5b08aa (feat(app): wire real web bridge with SSE progress transport)
=======
  return (await response.json()) as T;
>>>>>>> c2b9da5e (feat: add per-command result types and hooks across TypeScript layer)
}

export function createWebBridge(): TreeTimeBridge {
  return {
    version: () => getJson<VersionInfo>("version"),
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
=======
    datasets: () => getJson<DatasetInfo[]>("datasets"),
>>>>>>> b8625b9a (feat(app): wire datasets through bridge contract and implementations)
    ancestral: (args, onProgress) => postSse("ancestral", args, onProgress),
    clock: (args, onProgress) => postSse("clock", args, onProgress),
    timetree: (args, onProgress) => postSse("timetree", args, onProgress),
    mugration: (args, onProgress) => postSse("mugration", args, onProgress),
    optimize: (args, onProgress) => postSse("optimize", args, onProgress),
    prune: (args, onProgress) => postSse("prune", args, onProgress),
>>>>>>> 3e5b08aa (feat(app): wire real web bridge with SSE progress transport)
  };
}
