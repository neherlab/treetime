import { contextBridge, ipcRenderer } from "electron";
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

async function invokeCommand<T>(channel: string, args: unknown): Promise<T> {
  const result = await ipcRenderer.invoke(channel, JSON.stringify(args));
  return typeof result === "string" ? (JSON.parse(result) as T) : (result as T);
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
};

contextBridge.exposeInMainWorld("treetime", bridge);
