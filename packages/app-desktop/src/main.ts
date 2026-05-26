import { app, BrowserWindow, ipcMain, nativeTheme } from "electron";
import * as path from "path";

if (process.env.ELECTRON_DISABLE_SANDBOX === "1") {
  app.commandLine.appendSwitch("no-sandbox");
}

const projectRoot = process.env.TREETIME_PROJECT_ROOT;
if (projectRoot) {
  process.chdir(projectRoot);
}

function registerThemeHandler() {
  ipcMain.on("treetime:theme", (_event, theme: string) => {
    nativeTheme.themeSource = theme as "system" | "light" | "dark";
  });
}

function registerIpcHandlers() {
  // eslint-disable-next-line @typescript-eslint/no-require-imports
  const addon = require("@neherlab/app-napi");

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
          const parsed = JSON.parse(eventJson) as { type: string; data: unknown };
          event.sender.send(`treetime:${parsed.type}`, parsed.data);
        } catch {
          // Window closed during computation
        }
      });
    });
  }
}

function createWindow() {
  const win = new BrowserWindow({
    width: 1200,
    height: 800,
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      preload: path.join(__dirname, "preload.js"),
    },
  });

  if (process.env.VITE_DEV_SERVER_URL) {
    win.loadURL(process.env.VITE_DEV_SERVER_URL);
    win.webContents.openDevTools({ mode: "bottom" });
  } else {
    win.loadFile(path.join(__dirname, "../dist/index.html"));
  }
}

app.whenReady().then(() => {
  registerThemeHandler();
  registerIpcHandlers();
  createWindow();
});

app.on("window-all-closed", () => {
  app.quit();
});
