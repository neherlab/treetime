import * as path from "path";

import { parseBridgeEvent } from "@neherlab/app-contracts";
import * as addon from "@neherlab/app-napi";
import { app, BrowserWindow, ipcMain, nativeTheme } from "electron";

import { initDiagnostics } from "./diagnostics";

initDiagnostics("treetime-desktop");

if (process.env["ELECTRON_DISABLE_SANDBOX"] === "1") {
  app.commandLine.appendSwitch("no-sandbox");
}

const projectRoot = process.env["TREETIME_PROJECT_ROOT"];

if (projectRoot) {
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
  ipcMain.handle("treetime:version", () => {
    return addon.version();
  });

  ipcMain.handle("treetime:datasets", () => {
    return addon.datasets();
  });

  ipcMain.on("treetime:cancel", () => {
    addon.cancel();
  });

  const commands = ["ancestral", "clock", "timetree", "mugration", "optimize", "prune"] as const;

  for (const cmd of commands) {
    ipcMain.handle(`treetime:${cmd}`, (event: Electron.IpcMainInvokeEvent, argsJson: string) => {
      console.log(`[TreeTime IPC] ${cmd} called`);

      return addon[cmd](argsJson, (err: Error | null, eventJson: string) => {
        if (err || event.sender.isDestroyed()) return;

        try {
          const parsed = parseBridgeEvent(JSON.parse(eventJson));
          event.sender.send(`treetime:${parsed.type}`, parsed.data);
        } catch (error: unknown) {
          console.error(`[TreeTime IPC] ${cmd} event failed`, error);
        }
      });
    });
  }
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

  if (devServerUrl) {
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
