import { app, BrowserWindow, ipcMain } from "electron";
import * as path from "path";

if (process.env.ELECTRON_DISABLE_SANDBOX === "1") {
  app.commandLine.appendSwitch("no-sandbox");
}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
=======
const projectRoot = process.env.TREETIME_PROJECT_ROOT;
if (projectRoot) {
  process.chdir(projectRoot);
}

>>>>>>> d8153c8f (feat(desktop): Wire up progress events and cancellation in Electron IPC)
function registerIpcHandlers() {
  // eslint-disable-next-line @typescript-eslint/no-require-imports
  const addon = require("@neherlab/app-napi");

<<<<<<< HEAD
  ipcMain.handle("treetime:version", () => {
    return addon.version();
  });

  ipcMain.handle("treetime:datasets", () => {
    return addon.datasets();
  });

<<<<<<< HEAD
=======
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
=======
  ipcMain.on("treetime:cancel", () => {
    addon.cancel();
  });

>>>>>>> d8153c8f (feat(desktop): Wire up progress events and cancellation in Electron IPC)
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

<<<<<<< HEAD
>>>>>>> 9f689333 (fix(app-desktop): add missing version IPC handler and QueryProvider)
=======
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
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
  registerIpcHandlers();
  createWindow();
});

app.on("window-all-closed", () => {
  app.quit();
});
