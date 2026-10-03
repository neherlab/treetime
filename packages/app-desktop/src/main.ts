import { mkdirSync } from "node:fs";
import * as path from "node:path";
import { pathToFileURL } from "node:url";

import { zPickFilesRequest, zPickFolderRequest } from "@neherlab/app-contracts";
import { appPaths } from "@neherlab/app-napi";
import {
  app,
  BrowserWindow,
  dialog,
  ipcMain,
  MessageChannelMain,
  nativeTheme,
  net,
  protocol,
  shell,
  utilityProcess,
  webContents,
  type IpcMainEvent,
  type IpcMainInvokeEvent,
  type OpenDialogOptions,
  type UtilityProcess,
  type WebContents,
} from "electron";

import { APP_DIR_ENV, checkoutAppDir } from "./app-dir";
import { APP_SCHEME, APP_SCHEME_PRIVILEGES, APP_URL, resolveAppAsset } from "./app-scheme";
import { shouldRestart, stopReason } from "./backend-process";
import type { ControlReply, SaveRequest } from "./backend-protocol";
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
import { DIAGNOSTIC_DIR_ENV, initDiagnostics } from "./diagnostics";
import { confineNavigation, isTrustedSender } from "./security";
import type { SaveReply, SaveRunArchiveDialog, SaveRunFileDialog } from "./shell-protocol";

const projectRoot = process.env["TREETIME_PROJECT_ROOT"];

if (projectRoot !== undefined && projectRoot !== "") {
  process.chdir(projectRoot);
}

const devServerUrl = process.env["VITE_DEV_SERVER_URL"];

const devServer = devServerUrl !== undefined && devServerUrl !== "";

if (!app.isPackaged) {
  process.env[APP_DIR_ENV] = checkoutAppDir(process.env, devServer ? "dev" : "prod", process.cwd());
}

const paths = appPaths();

mkdirSync(paths.profileDir, { recursive: true });

app.setPath("userData", paths.profileDir);

app.setAppLogsPath(paths.logsDir);

const diagnosticDir = process.env[DIAGNOSTIC_DIR_ENV] ?? path.join(paths.logsDir, "diagnostics");

initDiagnostics("treetime-desktop", diagnosticDir);

if (process.env["ELECTRON_DISABLE_SANDBOX"] === "1") {
  app.commandLine.appendSwitch("no-sandbox");
}

const appUrl = devServer ? devServerUrl : APP_URL;

const appRoot = path.join(__dirname, "../dist");

protocol.registerSchemesAsPrivileged([{ scheme: APP_SCHEME, privileges: APP_SCHEME_PRIVILEGES }]);

async function main(): Promise<void> {
  await app.whenReady();
  protocol.handle(APP_SCHEME, serveAppAsset);
  const backend = new BackendProcess();
  registerIpcHandlers(backend);
  await createWindow();
}

class BackendProcess {
  private child: UtilityProcess;
  private readonly exitTimes: number[] = [];
  private readonly saves = new Map<number, (reply: ControlReply) => void>();
  private readonly restarted: Array<() => void> = [];
  private nextSave = 0;
  private failure: string | undefined;

  constructor() {
    this.child = this.spawn();
  }

  connect(contents: WebContents): void {
    if (this.failure !== undefined) {
      contents.send(BACKEND_STOPPED_CHANNEL, this.failure, false);

      return;
    }

    const channel = new MessageChannelMain();
    this.child.postMessage({ kind: "port" }, [channel.port1]);
    contents.postMessage(BACKEND_PORT_CHANNEL, null, [channel.port2]);
  }

  save(request: (seq: number) => SaveRequest): Promise<ControlReply> {
    const seq = this.nextSave;
    this.nextSave += 1;

    if (this.failure !== undefined) {
      return Promise.resolve({ kind: "error", seq, error: this.failure });
    }

    return new Promise((resolve) => {
      this.saves.set(seq, resolve);
      this.child.postMessage(request(seq), []);
    });
  }

  restart(): Promise<void> {
    const { promise, resolve } = Promise.withResolvers<void>();
    this.restarted.push(resolve);

    if (this.restarted.length > 1) {
      return promise;
    }

    if (this.failure === undefined) {
      this.child.kill();
    } else {
      this.failure = undefined;
      this.exitTimes.length = 0;
      this.respawn();
    }

    return promise;
  }

  private spawn(): UtilityProcess {
    const child = utilityProcess.fork(path.join(__dirname, "backend.js"), [], {
      serviceName: "TreeTime back end",
      cwd: process.cwd(),
      stdio: "inherit",
      env: { ...process.env, [DIAGNOSTIC_DIR_ENV]: diagnosticDir },
    });

    child.on("message", (reply: ControlReply) => {
      this.saves.get(reply.seq)?.(reply);
      this.saves.delete(reply.seq);
    });
    child.once("exit", (code) => {
      this.exited(code);
    });

    return child;
  }

  private exited(code: number): void {
    const requested = this.restarted.length > 0;
    const now = performance.now();

    if (!requested) {
      this.exitTimes.push(now);
    }

    const restart = requested || shouldRestart(this.exitTimes, now);
    const reason = stopReason(code, requested, restart);

    console.error(`[TreeTime] ${reason}`);

    for (const [seq, resolve] of this.saves) {
      resolve({ kind: "error", seq, error: reason });
    }

    this.saves.clear();

    for (const contents of appContents()) {
      contents.send(BACKEND_STOPPED_CHANNEL, reason, restart);
    }

    if (!restart) {
      this.failure = reason;

      return;
    }

    this.respawn();
  }

  private respawn(): void {
    setImmediate(() => {
      this.child = this.spawn();

      for (const contents of appContents()) {
        this.connect(contents);
      }

      for (const resolve of this.restarted.splice(0)) {
        resolve();
      }
    });
  }
}

function registerIpcHandlers(backend: BackendProcess): void {
  handle(PICK_FILES_CHANNEL, (event, request) => pickFiles(event, request));
  handle(PICK_FOLDER_CHANNEL, (event, request) => pickFolder(event, request));
  handle(RESTART_BACKEND_CHANNEL, async () => backend.restart());
  handle(SAVE_RUN_FILE_CHANNEL, async (event, { id, path, name }: SaveRunFileDialog) =>
    saveTo(event, name, (destination) =>
      backend.save((seq) => ({ kind: "save-file", seq, request: { id, path, destination } })),
    ),
  );
  handle(SAVE_RUN_ARCHIVE_CHANNEL, async (event, { id, name }: SaveRunArchiveDialog) =>
    saveTo(event, name, (destination) =>
      backend.save((seq) => ({ kind: "save-archive", seq, request: { id, destination } })),
    ),
  );
  listen(BACKEND_PORT_REQUEST_CHANNEL, (event) => {
    backend.connect(event.sender);
  });
  listen(THEME_CHANNEL, (_event, theme) => {
    if (isThemeSource(theme)) {
      nativeTheme.themeSource = theme;
    }
  });
}

// oxlint-disable-next-line typescript/no-unnecessary-type-parameters -- Args types the arguments of each handler, which Electron passes untyped
function handle<Args extends unknown[]>(
  channel: string,
  handler: (event: IpcMainInvokeEvent, ...args: Args) => unknown,
): void {
  ipcMain.handle(channel, (event, ...args: Args) => {
    if (!isTrustedSender(event.senderFrame, appUrl)) {
      throw new Error(`${channel} refused a message from a frame outside the application`);
    }

    return handler(event, ...args);
  });
}

function listen(channel: string, listener: (event: IpcMainEvent, ...args: unknown[]) => void): void {
  ipcMain.on(channel, (event, ...args: unknown[]) => {
    if (isTrustedSender(event.senderFrame, appUrl)) {
      listener(event, ...args);
    }
  });
}

async function serveAppAsset(request: Request): Promise<Response> {
  const asset = resolveAppAsset(request.url, appRoot);

  if (asset.kind === "not-found") {
    return new Response(null, { status: 404 });
  }

  try {
    return await net.fetch(pathToFileURL(asset.path).href, { signal: request.signal });
  } catch {
    return new Response(null, { status: 404 });
  }
}

function isThemeSource(value: unknown): value is "system" | "light" | "dark" {
  return value === "system" || value === "light" || value === "dark";
}

function appContents(): WebContents[] {
  return webContents.getAllWebContents().filter((contents) => isTrustedSender(contents.mainFrame, appUrl));
}

async function pickFiles(event: IpcMainInvokeEvent, data: unknown): Promise<string[]> {
  const request = zPickFilesRequest.parse(data);

  const options: OpenDialogOptions = {
    title: request.title,
    properties: request.multiple ? ["openFile", "multiSelections"] : ["openFile"],
    filters: request.extensions.length > 0 ? [{ name: request.title, extensions: request.extensions }] : [],
  };

  const window = BrowserWindow.fromWebContents(event.sender);
  const result = window === null ? await dialog.showOpenDialog(options) : await dialog.showOpenDialog(window, options);

  return result.canceled ? [] : result.filePaths;
}

async function pickFolder(event: IpcMainInvokeEvent, data: unknown): Promise<string | null> {
  const request = zPickFolderRequest.parse(data);
  const options: OpenDialogOptions = { title: request.title, properties: ["openDirectory", "createDirectory"] };

  const window = BrowserWindow.fromWebContents(event.sender);
  const result = window === null ? await dialog.showOpenDialog(options) : await dialog.showOpenDialog(window, options);

  return result.canceled ? null : (result.filePaths[0] ?? null);
}

async function saveTo(
  event: IpcMainInvokeEvent,
  name: string,
  save: (destination: string) => Promise<ControlReply>,
): Promise<SaveReply> {
  const window = BrowserWindow.fromWebContents(event.sender);
  const options = { defaultPath: name };
  const choice = window === null ? await dialog.showSaveDialog(options) : await dialog.showSaveDialog(window, options);

  if (choice.canceled || choice.filePath === "") {
    return { saved: false };
  }

  const reply = await save(choice.filePath);

  return reply.kind === "saved" ? { saved: true } : { error: reply.error };
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
