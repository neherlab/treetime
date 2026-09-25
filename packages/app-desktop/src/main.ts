import * as path from "path";

import { zCreateRunRequest, zPickFilesRequest, zRunRecord } from "@neherlab/app-contracts";
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

import { PICK_FILES_CHANNEL, RUN_EVENT_CHANNEL, type IpcReply } from "./desktop-bridge";
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
  const runs = new addon.RunService(path.join(app.getPath("userData"), "runs"));
  const subscriptions = new Map<string, addon.Subscription>();

  const startRun = (id: string, configJson: string | null) => {
    void runToEnd(runs, id, configJson);
  };

  handle(PICK_FILES_CHANNEL, (event, requestJson) => pickFiles(event, text(requestJson)));
  handle("treetime:version", () => addon.version());
  handle("treetime:datasets", () => addon.datasets());
  handle("treetime:check-config", (_event, requestJson) => addon.checkConfigJson(text(requestJson)));
  handle("treetime:run-config", (_event, requestJson) => addon.runConfigJson(text(requestJson)));
  handle("treetime:check-inputs", (_event, requestJson) => addon.checkInputsJson(text(requestJson)));
  handle("treetime:runs:list", () => runs.list());
  handle("treetime:runs:get", (_event, id) => runs.get(text(id)));
  handle("treetime:runs:create", (_event, requestJson) => {
    const request = text(requestJson);
    const record = zRunRecord.parse(JSON.parse(runs.create(request)));

    if (!zCreateRunRequest.parse(JSON.parse(request)).defer_start) {
      startRun(record.id, null);
    }

    return runs.get(record.id);
  });
  handle("treetime:runs:start", (_event, id, configJson) => {
    startRun(text(id), configJson === null ? null : text(configJson));

    return runs.get(text(id));
  });
  handle("treetime:runs:update", (_event, id, requestJson) => runs.update(text(id), text(requestJson)));
  handle("treetime:runs:cancel", (_event, id) => JSON.stringify({ cancelled: runs.cancel(text(id)) }));
  handle("treetime:runs:delete", (_event, id) => {
    runs.delete(text(id));
  });
  handle("treetime:runs:restore", (_event, id) => runs.restore(text(id)));
  handle("treetime:runs:purge", (_event, id) => {
    runs.purge(text(id));
  });
  handle("treetime:runs:files", (_event, id) => runs.files(text(id)));
  handle("treetime:runs:read-file", (_event, id, filePath) => runs.readFile(text(id), text(filePath)));
  handle("treetime:runs:archive", (_event, id) => runs.archive(text(id)));
  handle("treetime:runs:results", (_event, id) => runs.results(text(id)));
  handle("treetime:runs:compare", (_event, id, other) => runs.compare(text(id), text(other)));
  handle("treetime:runs:clade-in-runs", (_event, requestJson) => runs.cladeInRuns(text(requestJson)));
  handle("treetime:runs:subscribe", (event, subscriptionId, id, from) => {
    const key = text(subscriptionId);
    const sender = event.sender;

    const subscription = runs.subscribe(text(id), Number(from), (err: Error | null, eventJson: string) => {
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

async function runToEnd(runs: addon.RunService, id: string, configJson: string | null): Promise<void> {
  try {
    const terminalJson = await runs.start(id, configJson);
    console.log(`[TreeTime IPC] run ${id} ended: ${terminalJson}`);
  } catch (error: unknown) {
    console.error(`[TreeTime IPC] run ${id} failed to run`, error);
  }
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
