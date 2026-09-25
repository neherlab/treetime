import { contextBridge, ipcRenderer, webUtils, type IpcRendererEvent } from "electron";

import {
  BACKEND_PORT_CHANNEL,
  BACKEND_PORT_REQUEST_CHANNEL,
  BACKEND_STOPPED_CHANNEL,
  PICK_FILES_CHANNEL,
  THEME_CHANNEL,
} from "./channels";
import type { DesktopShell, ShellMessage } from "./desktop-bridge";

interface MainWorld {
  postMessage(message: ShellMessage, targetOrigin: string, transfer?: IpcRendererEvent["ports"]): void;
}

declare const window: MainWorld;

ipcRenderer.on(BACKEND_PORT_CHANNEL, (event) => {
  window.postMessage({ channel: BACKEND_PORT_CHANNEL }, "*", event.ports);
});

ipcRenderer.on(BACKEND_STOPPED_CHANNEL, (_event, reason: string) => {
  window.postMessage({ channel: BACKEND_STOPPED_CHANNEL, reason }, "*");
});

const shell: DesktopShell = {
  connectBackend: () => {
    ipcRenderer.send(BACKEND_PORT_REQUEST_CHANNEL);
  },
  pickFiles: (request) => ipcRenderer.invoke(PICK_FILES_CHANNEL, request),
  pathForFile: (file) => webUtils.getPathForFile(file),
};

contextBridge.exposeInMainWorld("treetimeShell", shell);

contextBridge.exposeInMainWorld("electronTheme", {
  setTheme: (theme: string) => {
    ipcRenderer.send(THEME_CHANNEL, theme);
  },
});
