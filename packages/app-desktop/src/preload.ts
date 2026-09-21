import { contextBridge, ipcRenderer } from "electron";

import { createDesktopBridge } from "./desktop-bridge";

contextBridge.exposeInMainWorld("treetime", createDesktopBridge(ipcRenderer));

contextBridge.exposeInMainWorld("electronTheme", {
  setTheme: (theme: string) => ipcRenderer.send("treetime:theme", theme),
});
