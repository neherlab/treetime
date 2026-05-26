import { contextBridge, ipcRenderer } from "electron";
import type { BridgeTransport, CommandOptions, LogEvent, ProgressEvent } from "@neherlab/app-contracts";
import { CancelledError, createBridge } from "@neherlab/app-contracts";

async function invokeQuery<T>(endpoint: string): Promise<T> {
  const json = await ipcRenderer.invoke(`treetime:${endpoint}`);
  return typeof json === "string" ? (JSON.parse(json) as T) : (json as T);
}

function logHandler(_event: Electron.IpcRendererEvent, data: LogEvent) {
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
}

function abortHandler() {
  ipcRenderer.send("treetime:cancel");
}

async function invokeCommand<T>(endpoint: string, args: unknown, options?: CommandOptions): Promise<T> {
  const progressHandler = (_event: Electron.IpcRendererEvent, data: ProgressEvent) => {
    options?.onProgress?.(data);
  };

  ipcRenderer.on("treetime:progress", progressHandler);
  ipcRenderer.on("treetime:log", logHandler);

  options?.signal?.addEventListener("abort", abortHandler);

  try {
    const result = await ipcRenderer.invoke(`treetime:${endpoint}`, JSON.stringify(args));
    return typeof result === "string" ? (JSON.parse(result) as T) : (result as T);
  } catch (err: unknown) {
    console.error("[TreeTime IPC]", endpoint, err);
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

function createDesktopTransport(): BridgeTransport {
  return {
    query: <T>(endpoint: string) => invokeQuery<T>(endpoint),
    command: <T>(endpoint: string, args: unknown, options?: CommandOptions) =>
      invokeCommand<T>(endpoint, args, options),
  };
}

contextBridge.exposeInMainWorld("treetime", createBridge(createDesktopTransport()));
contextBridge.exposeInMainWorld("electronTheme", {
  setTheme: (theme: string) => ipcRenderer.send("treetime:theme", theme),
});
