import { app, BrowserWindow, ipcMain } from "electron";
import * as path from "path";

if (process.env.ELECTRON_DISABLE_SANDBOX === "1") {
  app.commandLine.appendSwitch("no-sandbox");
}

<<<<<<< HEAD
=======
function registerIpcHandlers() {
  // eslint-disable-next-line @typescript-eslint/no-require-imports
  const addon = require("@neherlab/app-napi");

  ipcMain.handle("treetime:version", () => {
    return addon.version();
  });

  const commands = ["ancestral", "clock", "timetree", "mugration", "optimize", "prune"] as const;
  for (const cmd of commands) {
    ipcMain.handle(`treetime:${cmd}`, (_event: Electron.IpcMainInvokeEvent, argsJson: string) => {
      return addon[cmd](argsJson);
    });
  }
}

>>>>>>> 9f689333 (fix(app-desktop): add missing version IPC handler and QueryProvider)
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

  win.loadFile(path.join(__dirname, "..", "static", "index.html"));
}

app.whenReady().then(createWindow);

app.on("window-all-closed", () => {
  app.quit();
});
