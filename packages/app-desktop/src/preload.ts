import { contextBridge, ipcRenderer } from "electron";
import type { TreeTimeBridge } from "@neherlab/app-contracts";

const bridge: TreeTimeBridge = {
  ancestral: (args) => ipcRenderer.invoke("treetime:ancestral", JSON.stringify(args)),
  clock: (args) => ipcRenderer.invoke("treetime:clock", JSON.stringify(args)),
  timetree: (args) => ipcRenderer.invoke("treetime:timetree", JSON.stringify(args)),
  mugration: (args) => ipcRenderer.invoke("treetime:mugration", JSON.stringify(args)),
  optimize: (args) => ipcRenderer.invoke("treetime:optimize", JSON.stringify(args)),
  prune: (args) => ipcRenderer.invoke("treetime:prune", JSON.stringify(args)),
};

contextBridge.exposeInMainWorld("treetime", bridge);
