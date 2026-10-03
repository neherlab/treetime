import { contextBridge, ipcRenderer, webUtils, type IpcRendererEvent } from "electron";

import {
  BACKEND_PORT_CHANNEL,
  BACKEND_PORT_REQUEST_CHANNEL,
  BACKEND_STOPPED_CHANNEL,
  PICK_FILES_CHANNEL,
  PICK_FOLDER_CHANNEL,
  RESTART_BACKEND_CHANNEL,
  SAVE_RUN_ARCHIVE_CHANNEL,
  SAVE_RUN_FILE_CHANNEL,
  THEME_CHANNEL,
} from "./channels";
import type { DesktopShell } from "./desktop-shell";

interface MainWorld {
  postMessage(message: { channel: string }, targetOrigin: string, transfer?: IpcRendererEvent["ports"]): void;
}

declare const window: MainWorld;

ipcRenderer.on(BACKEND_PORT_CHANNEL, (event) => {
  window.postMessage({ channel: BACKEND_PORT_CHANNEL }, "*", event.ports);
});

const shell: DesktopShell = {
  connectBackend: () => {
    ipcRenderer.send(BACKEND_PORT_REQUEST_CHANNEL);
  },
  onBackendStopped: (listener) => {
    ipcRenderer.on(BACKEND_STOPPED_CHANNEL, (_event, reason: string, restarts: boolean) => {
      listener(reason, restarts);
    });
  },
  pickFiles: (request) => ipcRenderer.invoke(PICK_FILES_CHANNEL, request),
  pickFolder: (request) => ipcRenderer.invoke(PICK_FOLDER_CHANNEL, request),
  restartBackend: () => ipcRenderer.invoke(RESTART_BACKEND_CHANNEL),
  saveRunFile: (request) => ipcRenderer.invoke(SAVE_RUN_FILE_CHANNEL, request),
  saveRunArchive: (request) => ipcRenderer.invoke(SAVE_RUN_ARCHIVE_CHANNEL, request),
  pathForFile: (file) => webUtils.getPathForFile(file),
};

contextBridge.exposeInMainWorld("treetimeShell", shell);

contextBridge.exposeInMainWorld("electronTheme", {
  setTheme: (theme: string) => {
    ipcRenderer.send(THEME_CHANNEL, theme);
  },
});
