import * as path from "path";

import { parseJobEvent } from "@neherlab/app-contracts";
import * as addon from "@neherlab/app-napi";
import { app, BrowserWindow, ipcMain, nativeTheme } from "electron";

import { JOB_EVENT_CHANNEL } from "./desktop-bridge";
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
  const runner = new addon.CommandRunner();

  ipcMain.handle("treetime:version", () => {
    return addon.version();
  });

  ipcMain.handle("treetime:datasets", () => {
    return addon.datasets();
  });

  ipcMain.handle("treetime:check-config", (_event: Electron.IpcMainInvokeEvent, requestJson: string) => {
    return addon.checkConfigJson(requestJson);
  });

  ipcMain.on("treetime:cancel", (_event: Electron.IpcMainEvent, jobId: string) => {
    runner.cancel(jobId);
  });

  ipcMain.handle(
    "treetime:run",
    (event: Electron.IpcMainInvokeEvent, jobId: string, command: string, configJson: string) => {
      console.log(`[TreeTime IPC] ${command} job ${jobId} started`);

      return runner.run(jobId, command, configJson, (err: Error | null, eventJson: string) => {
        if (err || event.sender.isDestroyed()) return;

        try {
          parseJobEvent(JSON.parse(eventJson));
          event.sender.send(JOB_EVENT_CHANNEL, jobId, eventJson);
        } catch (error: unknown) {
          console.error(`[TreeTime IPC] ${command} job ${jobId} event failed`, error);
        }
      });
    },
  );
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
