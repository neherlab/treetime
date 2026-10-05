import { mkdirSync } from "node:fs";
import * as path from "node:path";
import { pathToFileURL } from "node:url";

import { errorMessage, zUiTheme, type ErrorResponse, type UiTheme } from "@neherlab/app-contracts";
import { ApiError, createApiClient, runsSave, type ApiClient } from "@neherlab/app-contracts/client";
import { appStartup } from "@neherlab/app-napi";
import type { BackendStopped, HostHandlerReply, HostHandlerRequest } from "@neherlab/app-ui/host";
import { windowBackground, WINDOW_MIN_SIZE } from "@neherlab/app-ui/window";
import {
  app,
  BrowserWindow,
  dialog,
  ipcMain,
  MessageChannelMain,
  nativeTheme,
  net,
  protocol,
  screen,
  shell,
  utilityProcess,
  webContents,
  type MessagePortMain,
  type OpenDialogOptions,
  type OpenDialogReturnValue,
  type UtilityProcess,
  type WebContents,
} from "electron";

import { appFolderEnv } from "./app-dir";
import { APP_SCHEME, APP_SCHEME_PRIVILEGES, APP_URL, resolveAppAsset } from "./app-scheme";
import { backendStop } from "./backend-process";
import type { ControlReply, ControlRequest, PortScope } from "./backend-protocol";
import { DIAGNOSTIC_DIR_ENV, initDiagnostics } from "./diagnostics";
import { HostConnection, mainFetchPort } from "./host-port";
import { emit, handle, listenForPortRequests, sendPort } from "./ipc-main";
import { napiErrorResponse } from "./napi-error";
import { createPortFetch } from "./port-fetch";
import { confineNavigation, isTrustedSender } from "./security";
import { firstWindowBounds } from "./window-bounds";

const MAIN_WINDOW_NAME = "main";

const HOST_BASE_URL = "http://treetime.host";

const launchDir = process.cwd();

const appDirSwitch = app.commandLine.getSwitchValue("app-dir");

const projectRoot = process.env["TREETIME_PROJECT_ROOT"];

if (projectRoot !== undefined && projectRoot !== "") {
  process.chdir(projectRoot);
}

const devServerUrl = process.env["VITE_DEV_SERVER_URL"];

const devServer = devServerUrl !== undefined && devServerUrl !== "";

Object.assign(
  process.env,
  appFolderEnv({
    appDirSwitch,
    env: process.env,
    launchDir,
    checkout: app.isPackaged ? undefined : { mode: devServer ? "dev" : "prod", root: process.cwd() },
  }),
);

if (process.env["ELECTRON_DISABLE_SANDBOX"] === "1") {
  app.commandLine.appendSwitch("no-sandbox");
}

const appUrl = devServer ? devServerUrl : APP_URL;

const appRoot = path.join(__dirname, "../dist");

protocol.registerSchemesAsPrivileged([{ scheme: APP_SCHEME, privileges: APP_SCHEME_PRIVILEGES }]);

function readStartup(): Startup | undefined {
  try {
    const { profileDir, logsDir, theme } = appStartup();

    return { profileDir, logsDir, theme: zUiTheme.parse(theme) };
  } catch (error: unknown) {
    const { message, causes } = napiErrorResponse(error);
    dialog.showErrorBox("TreeTime cannot start", [message, ...causes].join(": "));
    app.quit();

    return undefined;
  }
}

function launch({ logsDir, theme }: Startup): void {
  const diagnosticDir = process.env[DIAGNOSTIC_DIR_ENV] ?? path.join(logsDir, "diagnostics");
  initDiagnostics("treetime-desktop", diagnosticDir);
  const mainWindow = new MainWindow();

  app.on("second-instance", () => {
    void mainWindow.show();
  });
  app.on("window-all-closed", () => {
    if (process.platform !== "darwin") {
      app.quit();
    }
  });
  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- a failed start rejects with an untyped error, logged before the app quits
  main(mainWindow, theme, diagnosticDir).catch((error: unknown) => {
    console.error(error);
    app.quit();
  });
}

async function main(mainWindow: MainWindow, theme: UiTheme, diagnosticDir: string): Promise<void> {
  await app.whenReady();
  protocol.handle(APP_SCHEME, serveAppAsset);
  nativeTheme.themeSource = theme;
  nativeTheme.on("updated", () => {
    for (const window of BrowserWindow.getAllWindows()) {
      window.setBackgroundColor(windowBackground(nativeTheme.shouldUseDarkColors));
    }
  });
  const backend = new BackendProcess(diagnosticDir);
  registerIpcHandlers(backend, createApiClient({ baseUrl: HOST_BASE_URL, fetch: createPortFetch(backend.host) }));
  await mainWindow.show();
  app.on("activate", () => {
    void mainWindow.show();
  });
}

interface Startup {
  profileDir: string;
  logsDir: string;
  theme: UiTheme;
}

class MainWindow {
  private window: BrowserWindow | undefined;

  async show(): Promise<void> {
    await app.whenReady();

    if (this.window === undefined || this.window.isDestroyed()) {
      this.window = createWindow();

      return;
    }

    if (this.window.isMinimized()) {
      this.window.restore();
    }

    this.window.show();
    this.window.focus();
  }
}

class BackendProcess {
  readonly host = new HostConnection();
  private child: UtilityProcess;
  private readonly exitTimes: number[] = [];
  private readonly restarted: Array<() => void> = [];
  private failure: BackendStopped | undefined;
  private startError: ErrorResponse | undefined;
  private readonly diagnosticDir: string;

  constructor(diagnosticDir: string) {
    this.diagnosticDir = diagnosticDir;
    this.child = this.spawn();
    this.connectHost();
  }

  connect(contents: WebContents): void {
    if (this.failure !== undefined) {
      emit(contents, "backend-stopped", this.failure);

      return;
    }

    sendPort(contents, this.openPort("renderer"));
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

  private connectHost(): void {
    this.host.connect(mainFetchPort(this.openPort("host")));
  }

  private openPort(scope: PortScope): MessagePortMain {
    const channel = new MessageChannelMain();
    this.child.postMessage({ kind: "port", scope } satisfies ControlRequest, [channel.port1]);

    return channel.port2;
  }

  private spawn(): UtilityProcess {
    const child = utilityProcess.fork(path.join(__dirname, "backend.js"), [], {
      serviceName: "TreeTime back end",
      cwd: process.cwd(),
      stdio: "inherit",
      env: { ...process.env, [DIAGNOSTIC_DIR_ENV]: this.diagnosticDir },
    });

    this.startError = undefined;
    child.on("message", (reply: ControlReply) => {
      this.startError = reply.error;
      child.kill();
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

    const stop = backendStop({ code, requested, startError: this.startError, exitTimes: this.exitTimes, now });

    console.error(`[TreeTime] ${stop.reason}`);
    this.host.stop(stop);

    for (const contents of appContents()) {
      emit(contents, "backend-stopped", stop);
    }

    if (!stop.restarts) {
      this.failure = stop;

      return;
    }

    this.respawn();
  }

  private respawn(): void {
    setImmediate(() => {
      this.child = this.spawn();
      this.connectHost();

      for (const contents of appContents()) {
        this.connect(contents);
      }

      for (const resolve of this.restarted.splice(0)) {
        resolve();
      }
    });
  }
}

function registerIpcHandlers(backend: BackendProcess, hostClient: ApiClient): void {
  handle(ipcMain, appUrl, "pick-files", async (sender, request) => pickFiles(sender, request));
  handle(ipcMain, appUrl, "pick-folder", async (sender, request) => pickFolder(sender, request));
  handle(ipcMain, appUrl, "save-run", async (sender, request) => saveRun(sender, request, hostClient));
  handle(ipcMain, appUrl, "restart-backend", async () => {
    await backend.restart();

    return undefined;
  });
  handle(ipcMain, appUrl, "set-native-theme", (_sender, theme) => {
    nativeTheme.themeSource = theme;

    return undefined;
  });
  listenForPortRequests(ipcMain, appUrl, (sender) => {
    backend.connect(sender);
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

function appContents(): WebContents[] {
  return webContents.getAllWebContents().filter((contents) => isTrustedSender(contents.mainFrame, appUrl));
}

async function pickFiles(
  sender: WebContents,
  request: HostHandlerRequest<"pick-files">,
): Promise<HostHandlerReply<"pick-files">> {
  const options: OpenDialogOptions = {
    title: request.title,
    properties: request.multiple ? ["openFile", "multiSelections"] : ["openFile"],
    filters: request.extensions.length > 0 ? [{ name: request.title, extensions: request.extensions }] : [],
  };

  const result = await showOpenDialog(sender, options);

  return result.canceled ? [] : result.filePaths;
}

async function pickFolder(
  sender: WebContents,
  request: HostHandlerRequest<"pick-folder">,
): Promise<HostHandlerReply<"pick-folder">> {
  const result = await showOpenDialog(sender, {
    title: request.title,
    properties: ["openDirectory", "createDirectory"],
  });

  return result.canceled ? undefined : result.filePaths[0];
}

async function showOpenDialog(sender: WebContents, options: OpenDialogOptions): Promise<OpenDialogReturnValue> {
  const window = BrowserWindow.fromWebContents(sender);

  return window === null ? dialog.showOpenDialog(options) : dialog.showOpenDialog(window, options);
}

async function saveRun(
  sender: WebContents,
  { id, path: file, name }: HostHandlerRequest<"save-run">,
  client: ApiClient,
): Promise<HostHandlerReply<"save-run">> {
  const window = BrowserWindow.fromWebContents(sender);
  const options = { defaultPath: name };
  const choice = window === null ? await dialog.showSaveDialog(options) : await dialog.showSaveDialog(window, options);

  if (choice.canceled || choice.filePath === "") {
    return { kind: "canceled" };
  }

  try {
    await runsSave({ client, path: { id }, body: { path: file, destination: choice.filePath }, throwOnError: true });

    return { kind: "saved" };
  } catch (error: unknown) {
    return { kind: "error", error: errorResponse(error) };
  }
}

// oxlint-disable-next-line anti-slop/no-unknown-parameters -- a rejected request throws an ApiError or a transport error, sorted here
function errorResponse(error: unknown): ErrorResponse {
  return error instanceof ApiError
    ? error.response
    : { code: "internal_error", message: errorMessage(error), causes: [] };
}

function createWindow(): BrowserWindow {
  const bounds = firstWindowBounds(screen.getPrimaryDisplay().workArea, WINDOW_MIN_SIZE);

  const win = new BrowserWindow({
    ...bounds,
    name: MAIN_WINDOW_NAME,
    windowStatePersistence: true,
    minWidth: WINDOW_MIN_SIZE.width,
    minHeight: WINDOW_MIN_SIZE.height,
    show: false,
    backgroundColor: windowBackground(nativeTheme.shouldUseDarkColors),
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      sandbox: true,
      preload: path.join(__dirname, "preload.js"),
    },
  });

  win.once("ready-to-show", () => {
    win.show();
  });
  confineNavigation(win.webContents, appUrl, (url) => void shell.openExternal(url));
  void loadWindow(win);

  return win;
}

async function loadWindow(win: BrowserWindow): Promise<void> {
  try {
    await win.loadURL(appUrl);
  } catch (error: unknown) {
    console.error(error);

    return;
  }

  if (devServer) {
    win.webContents.openDevTools();
  }
}

const startup = readStartup();

if (startup !== undefined) {
  mkdirSync(startup.profileDir, { recursive: true });
  app.setPath("userData", startup.profileDir);
  app.setAppLogsPath(startup.logsDir);

  if (app.requestSingleInstanceLock()) {
    launch(startup);
  } else {
    app.quit();
  }
}
