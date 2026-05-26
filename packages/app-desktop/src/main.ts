import { app, BrowserWindow, ipcMain } from "electron";
import * as path from "path";

if (process.env.ELECTRON_DISABLE_SANDBOX === "1") {
  app.commandLine.appendSwitch("no-sandbox");
}

<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
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

=======
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
  const commands = ["ancestral", "clock", "timetree", "mugration", "optimize", "prune"] as const;
  for (const cmd of commands) {
    ipcMain.handle(`treetime:${cmd}`, (_event: Electron.IpcMainInvokeEvent, argsJson: string) => {
      return addon[cmd](argsJson);
    });
  }
}

<<<<<<< HEAD
>>>>>>> 9f689333 (fix(app-desktop): add missing version IPC handler and QueryProvider)
=======
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
function createWindow() {
  const win = new BrowserWindow({
    width: 900,
    height: 600,
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      preload: path.join(__dirname, "preload.js"),
    },
  });

  const viteDevUrl = process.env.VITE_DEV_SERVER_URL;
  if (viteDevUrl) {
    win.loadURL(viteDevUrl);
  } else {
    win.loadFile(path.join(__dirname, "..", "renderer", "index.html"));
  }
}

app.whenReady().then(() => {
  registerIpcHandlers();
  createWindow();
});

app.on("window-all-closed", () => {
  app.quit();
});
