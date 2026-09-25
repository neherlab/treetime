import * as path from "path";

import { parseRunEvent, zCreateRunRequest, zRunRecord } from "@neherlab/app-contracts";
import * as addon from "@neherlab/app-napi";
import { app, BrowserWindow, ipcMain, nativeTheme } from "electron";

import { RUN_EVENT_CHANNEL } from "./desktop-bridge";
import { initDiagnostics } from "./diagnostics";

initDiagnostics("treetime-desktop");

if (process.env["ELECTRON_DISABLE_SANDBOX"] === "1") {
  app.commandLine.appendSwitch("no-sandbox");
}

const projectRoot = process.env["TREETIME_PROJECT_ROOT"];

if (projectRoot !== undefined && projectRoot !== "") {
  process.chdir(projectRoot);
}

async function main(): Promise<void> {
  await app.whenReady();
  registerThemeHandler();
  registerIpcHandlers();
  await createWindow();
}

function registerThemeHandler(): void {
  ipcMain.on("treetime:theme", (_event, theme: string) => {
    if (isThemeSource(theme)) {
      nativeTheme.themeSource = theme;
    }
  });
}

function isThemeSource(value: string): value is "system" | "light" | "dark" {
  return value === "system" || value === "light" || value === "dark";
}

function registerIpcHandlers(): void {
  const runs = new addon.RunService(path.join(app.getPath("userData"), "runs"));

  const startRun = (id: string, configJson: string | null) => {
    void runToEnd(runs, id, configJson);
  };

  ipcMain.handle("treetime:version", () => addon.version());
  ipcMain.handle("treetime:datasets", () => addon.datasets());
  ipcMain.handle("treetime:check-config", (_event, requestJson: string) => addon.checkConfigJson(requestJson));
  ipcMain.handle("treetime:run-config", (_event, requestJson: string) => addon.runConfigJson(requestJson));
  ipcMain.handle("treetime:check-inputs", (_event, requestJson: string) => addon.checkInputsJson(requestJson));
  ipcMain.handle("treetime:runs:list", () => runs.list());
  ipcMain.handle("treetime:runs:get", (_event, id: string) => runs.get(id));
  ipcMain.handle("treetime:runs:create", (_event, requestJson: string) => {
    const record = parseRecord(runs.create(requestJson));

    if (!parseDeferStart(requestJson)) {
      startRun(record.id, null);
    }

    return runs.get(record.id);
  });
  ipcMain.handle("treetime:runs:start", (_event, id: string, configJson: string | null) => {
    startRun(id, configJson);

    return runs.get(id);
  });
  ipcMain.handle("treetime:runs:update", (_event, id: string, requestJson: string) => runs.update(id, requestJson));
  ipcMain.handle("treetime:runs:cancel", (_event, id: string) => JSON.stringify({ cancelled: runs.cancel(id) }));
  ipcMain.handle("treetime:runs:delete", (_event, id: string) => {
    runs.delete(id);
  });
  ipcMain.handle("treetime:runs:restore", (_event, id: string) => runs.restore(id));
  ipcMain.handle("treetime:runs:purge", (_event, id: string) => {
    runs.purge(id);
  });
  ipcMain.handle("treetime:runs:files", (_event, id: string) => runs.files(id));
  ipcMain.handle("treetime:runs:read-file", (_event, id: string, filePath: string) => runs.readFile(id, filePath));
  ipcMain.handle("treetime:runs:archive", (_event, id: string) => runs.archive(id));
  ipcMain.handle(
    "treetime:runs:subscribe",
    (event: Electron.IpcMainInvokeEvent, subscriptionId: string, id: string, from: number) => {
      runs.subscribe(id, from, (err: Error | null, eventJson: string) => {
        if (err || event.sender.isDestroyed()) return;

        try {
          parseRunEvent(JSON.parse(eventJson));
          event.sender.send(RUN_EVENT_CHANNEL, subscriptionId, eventJson);
        } catch (error: unknown) {
          console.error(`[TreeTime IPC] event of run ${id} is malformed`, error);
        }
      });
    },
  );
}

async function runToEnd(runs: addon.RunService, id: string, configJson: string | null): Promise<void> {
  try {
    const terminalJson = await runs.start(id, configJson);
    console.log(`[TreeTime IPC] run ${id} ended: ${terminalJson}`);
  } catch (error: unknown) {
    console.error(`[TreeTime IPC] run ${id} failed to run`, error);
  }
}

function parseRecord(json: string): { id: string } {
  return zRunRecord.parse(JSON.parse(json));
}

function parseDeferStart(requestJson: string): boolean {
  return zCreateRunRequest.parse(JSON.parse(requestJson)).defer_start;
}

async function createWindow(): Promise<void> {
  const win = new BrowserWindow({
    width: 1200,
    height: 800,
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      preload: path.join(__dirname, "preload.js"),
    },
  });

  const devServerUrl = process.env["VITE_DEV_SERVER_URL"];

  if (devServerUrl !== undefined && devServerUrl !== "") {
    await win.loadURL(devServerUrl);
    win.webContents.openDevTools({ mode: "bottom" });
  } else {
    await win.loadFile(path.join(__dirname, "../dist/index.html"));
  }
}

main().catch((error: unknown) => {
  console.error(error);
  app.quit();
});

app.on("window-all-closed", () => {
  app.quit();
});
