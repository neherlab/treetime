import { contextBridge, ipcRenderer, webUtils, type IpcRendererEvent } from "electron";

import { createPreloadHost, relayBackendPort, type PortTarget } from "./ipc-renderer";

declare const window: PortTarget<IpcRendererEvent["ports"][number]>;

relayBackendPort(ipcRenderer, window);

contextBridge.exposeInMainWorld("treetimeHost", createPreloadHost(ipcRenderer, webUtils));
