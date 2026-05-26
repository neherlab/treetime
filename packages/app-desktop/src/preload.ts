import { contextBridge, ipcRenderer } from "electron";
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
import type { TreeTimeBridge } from "@neherlab/app-contracts";

const bridge: TreeTimeBridge = {
  ancestral: (args) => ipcRenderer.invoke("treetime:ancestral", JSON.stringify(args)),
  clock: (args) => ipcRenderer.invoke("treetime:clock", JSON.stringify(args)),
  timetree: (args) => ipcRenderer.invoke("treetime:timetree", JSON.stringify(args)),
  mugration: (args) => ipcRenderer.invoke("treetime:mugration", JSON.stringify(args)),
  optimize: (args) => ipcRenderer.invoke("treetime:optimize", JSON.stringify(args)),
  prune: (args) => ipcRenderer.invoke("treetime:prune", JSON.stringify(args)),
=======
import type { TreeTimeBridge, CommandResult, VersionInfo } from "@neherlab/app-contracts";
=======
import type { TreeTimeBridge, VersionInfo } from "@neherlab/app-contracts";
>>>>>>> c2b9da5e (feat: add per-command result types and hooks across TypeScript layer)
=======
import type { TreeTimeBridge, DatasetInfo, VersionInfo } from "@neherlab/app-contracts";
>>>>>>> b8625b9a (feat(app): wire datasets through bridge contract and implementations)
=======
import type {
  CommandOptions,
  DatasetInfo,
  LogEvent,
  ProgressEvent,
  TreeTimeBridge,
  VersionInfo,
} from "@neherlab/app-contracts";
>>>>>>> d8153c8f (feat(desktop): Wire up progress events and cancellation in Electron IPC)

class CancelledError extends Error {
  constructor() {
    super("Operation cancelled");
    this.name = "CancelledError";
  }
}

async function invokeCommand<T>(channel: string, args: unknown, options?: CommandOptions): Promise<T> {
  const progressHandler = (_event: Electron.IpcRendererEvent, data: ProgressEvent) => {
    options?.onProgress?.(data);
  };

  const logHandler = (_event: Electron.IpcRendererEvent, data: LogEvent) => {
    switch (data.level) {
      case "Error":
        console.error(`[TreeTime] ${data.message}`);
        break;
      case "Warn":
        console.warn(`[TreeTime] ${data.message}`);
        break;
      default:
        console.log(`[TreeTime] [${data.level}] ${data.message}`);
        break;
    }
  };

  ipcRenderer.on("treetime:progress", progressHandler);
  ipcRenderer.on("treetime:log", logHandler);

  const abortHandler = () => {
    ipcRenderer.send("treetime:cancel");
  };
  options?.signal?.addEventListener("abort", abortHandler);

  try {
    const result = await ipcRenderer.invoke(channel, JSON.stringify(args));
    return typeof result === "string" ? (JSON.parse(result) as T) : (result as T);
  } catch (err: unknown) {
    if (err instanceof Error && err.message.includes("cancelled")) {
      throw new CancelledError();
    }
    throw err;
  } finally {
    ipcRenderer.removeListener("treetime:progress", progressHandler);
    ipcRenderer.removeListener("treetime:log", logHandler);
    options?.signal?.removeEventListener("abort", abortHandler);
  }
}

const bridge: TreeTimeBridge = {
  version: async (): Promise<VersionInfo> => {
    const json = await ipcRenderer.invoke("treetime:version");
    return JSON.parse(json) as VersionInfo;
  },
  datasets: async (): Promise<DatasetInfo[]> => {
    const json = await ipcRenderer.invoke("treetime:datasets");
    return JSON.parse(json) as DatasetInfo[];
  },
<<<<<<< HEAD
=======
import type { TreeTimeBridge, CommandResult } from "@neherlab/app-contracts";

async function invokeCommand(channel: string, args: unknown): Promise<CommandResult> {
  try {
    await ipcRenderer.invoke(channel, JSON.stringify(args));
    return { status: "ok" };
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : String(err);
    return { status: "error", error: message };
  }
}

const bridge: TreeTimeBridge = {
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
  ancestral: (args) => invokeCommand("treetime:ancestral", args),
  clock: (args) => invokeCommand("treetime:clock", args),
  timetree: (args) => invokeCommand("treetime:timetree", args),
  mugration: (args) => invokeCommand("treetime:mugration", args),
  optimize: (args) => invokeCommand("treetime:optimize", args),
  prune: (args) => invokeCommand("treetime:prune", args),
<<<<<<< HEAD
>>>>>>> 33bee034 (feat: add end-to-end version info across all layers)
=======
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
=======
  ancestral: (args, options) => invokeCommand("treetime:ancestral", args, options),
  clock: (args, options) => invokeCommand("treetime:clock", args, options),
  timetree: (args, options) => invokeCommand("treetime:timetree", args, options),
  mugration: (args, options) => invokeCommand("treetime:mugration", args, options),
  optimize: (args, options) => invokeCommand("treetime:optimize", args, options),
  prune: (args, options) => invokeCommand("treetime:prune", args, options),
>>>>>>> d8153c8f (feat(desktop): Wire up progress events and cancellation in Electron IPC)
};

contextBridge.exposeInMainWorld("treetime", bridge);
