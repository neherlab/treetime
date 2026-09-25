import * as path from "path";

import { zPickFilesRequest } from "@neherlab/app-contracts";
import * as addon from "@neherlab/app-napi";
import {
  app,
  BrowserWindow,
  dialog,
  ipcMain,
  nativeTheme,
  shell,
  type IpcMainEvent,
  type IpcMainInvokeEvent,
  type OpenDialogOptions,
} from "electron";

import { CALL_CHANNEL, PICK_FILES_CHANNEL, RUN_EVENT_CHANNEL, type IpcReply } from "./desktop-bridge";
import { initDiagnostics } from "./diagnostics";
import { confineNavigation, isTrustedSender } from "./security";

initDiagnostics("treetime-desktop");

if (process.env["ELECTRON_DISABLE_SANDBOX"] === "1") {
  app.commandLine.appendSwitch("no-sandbox");
}

const projectRoot = process.env["TREETIME_PROJECT_ROOT"];

if (projectRoot !== undefined && projectRoot !== "") {
  process.chdir(projectRoot);
}

const devServerUrl = process.env["VITE_DEV_SERVER_URL"];

const appUrl =
  devServerUrl !== undefined && devServerUrl !== ""
    ? devServerUrl
    : new URL(`file://${path.join(__dirname, "../dist/index.html")}`).href;

async function main(): Promise<void> {
  await app.whenReady();
  registerThemeHandler();
  registerIpcHandlers();
  await createWindow();
}

function registerThemeHandler(): void {
  listen("treetime:theme", (_event, theme) => {
    if (isThemeSource(theme)) {
      nativeTheme.themeSource = theme;
    }
  });
}

function isThemeSource(value: unknown): value is "system" | "light" | "dark" {
  return value === "system" || value === "light" || value === "dark";
}

function registerIpcHandlers(): void {
  const backend = new addon.Backend(path.join(app.getPath("userData"), "runs"));
  const subscriptions = new Map<string, addon.Subscription>();

  handle(PICK_FILES_CHANNEL, (event, requestJson) => pickFiles(event, text(requestJson)));
  handle(CALL_CHANNEL, (_event, requestJson) => backend.call(text(requestJson)));
  handle("treetime:runs:read-file", (_event, id, filePath) => backend.readFile(text(id), text(filePath)));
  handle("treetime:runs:archive", (_event, id) => backend.archive(text(id)));
  handle("treetime:runs:subscribe", (event, subscriptionId, id, from) => {
    const key = text(subscriptionId);
    const sender = event.sender;

    const subscription = backend.subscribe(text(id), Number(from), (err: Error | null, eventJson: string) => {
      if (err === null && !sender.isDestroyed()) {
        sender.send(RUN_EVENT_CHANNEL, key, eventJson);
      }
    });

    subscriptions.set(key, subscription);
    sender.once("destroyed", () => {
      subscription.unsubscribe();
      subscriptions.delete(key);
    });
  });
  listen("treetime:runs:unsubscribe", (_event, subscriptionId) => {
    const key = text(subscriptionId);
    subscriptions.get(key)?.unsubscribe();
    subscriptions.delete(key);
  });
}

function handle(channel: string, handler: (event: IpcMainInvokeEvent, ...args: unknown[]) => unknown): void {
  ipcMain.handle(channel, async (event, ...args: unknown[]): Promise<IpcReply> => {
    if (!isTrustedSender(event.senderFrame, appUrl)) {
      throw new Error(`${channel} refused a message from a frame outside the application`);
    }

    try {
      return { ok: true, value: await handler(event, ...args) };
    } catch (error: unknown) {
      return { ok: false, error: error instanceof Error ? error.message : String(error) };
    }
  });
}

function listen(channel: string, listener: (event: IpcMainEvent, ...args: unknown[]) => void): void {
  ipcMain.on(channel, (event, ...args: unknown[]) => {
    if (isTrustedSender(event.senderFrame, appUrl)) {
      listener(event, ...args);
    }
  });
}

function text(value: unknown): string {
  if (typeof value !== "string") {
    throw new TypeError(`expected a string argument, got ${typeof value}`);
  }

  return value;
}

async function pickFiles(event: IpcMainInvokeEvent, requestJson: string): Promise<string[]> {
  const request = zPickFilesRequest.parse(JSON.parse(requestJson));

  const options: OpenDialogOptions = {
    title: request.title,
    properties: request.multiple ? ["openFile", "multiSelections"] : ["openFile"],
    filters: request.extensions.length > 0 ? [{ name: request.title, extensions: request.extensions }] : [],
  };

  const window = BrowserWindow.fromWebContents(event.sender);
  const result = window === null ? await dialog.showOpenDialog(options) : await dialog.showOpenDialog(window, options);

  return result.canceled ? [] : result.filePaths;
}

async function createWindow(): Promise<void> {
  const win = new BrowserWindow({
    width: 1200,
    height: 800,
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      sandbox: true,
      preload: path.join(__dirname, "preload.js"),
    },
  });

  confineNavigation(win.webContents, appUrl, (url) => void shell.openExternal(url));

  if (devServerUrl !== undefined && devServerUrl !== "") {
    await win.loadURL(devServerUrl);
    win.webContents.openDevTools({ mode: "bottom" });
  } else {
    await win.loadURL(appUrl);
  }
}

main().catch((error: unknown) => {
  console.error(error);
  app.quit();
});

app.on("window-all-closed", () => {
  app.quit();
});
