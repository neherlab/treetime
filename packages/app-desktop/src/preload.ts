import { contextBridge, ipcRenderer, webUtils } from "electron";

import { createDesktopBridge, createLocalFiles } from "./desktop-bridge";

contextBridge.exposeInMainWorld("treetime", createDesktopBridge(ipcRenderer));

contextBridge.exposeInMainWorld(
  "treetimeFiles",
  createLocalFiles(ipcRenderer, (file) => webUtils.getPathForFile(file)),
);

contextBridge.exposeInMainWorld("electronTheme", {
  setTheme: (theme: string) => ipcRenderer.send("treetime:theme", theme),
});
