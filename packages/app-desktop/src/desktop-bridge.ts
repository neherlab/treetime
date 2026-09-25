import {
  createBridge,
  parseRunEvent,
  zPickedFiles,
  type BridgeTransport,
  type LocalFiles,
  type LogEvent,
  type TransportEventOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";

export interface IpcRendererLike {
  invoke(channel: string, ...args: unknown[]): Promise<unknown>;
  on(channel: string, listener: (event: unknown, ...args: unknown[]) => void): void;
  removeListener(channel: string, listener: (event: unknown, ...args: unknown[]) => void): void;
  send(channel: string, ...args: unknown[]): void;
}

export const RUN_EVENT_CHANNEL = "treetime:run-event";

export const PICK_FILES_CHANNEL = "treetime:pick-files";

export class LocalInputsError extends Error {
  constructor() {
    super("the desktop application reads inputs from local file paths; name the files in the run configuration");
    this.name = "LocalInputsError";
  }
}

export function createDesktopBridge(
  ipc: IpcRendererLike,
  newSubscriptionId: () => string = randomSubscriptionId,
): TreeTimeBridge {
  return createBridge(createDesktopTransport(ipc, newSubscriptionId));
}

export function createLocalFiles(ipc: IpcRendererLike, pathForFile: (file: File) => string): LocalFiles {
  return {
    async pickFiles(request) {
      return zPickedFiles.parse(await ipc.invoke(PICK_FILES_CHANNEL, JSON.stringify(request)));
    },
    pathForFile,
  };
}

function createDesktopTransport(ipc: IpcRendererLike, newSubscriptionId: () => string): BridgeTransport {
  async function call(channel: string, ...args: unknown[]): Promise<unknown> {
    return decode(await ipc.invoke(`treetime:${channel}`, ...args));
  }

  async function bytes(channel: string, ...args: unknown[]): Promise<Uint8Array> {
    const value = await ipc.invoke(`treetime:${channel}`, ...args);

    if (!(value instanceof Uint8Array)) {
      throw new TypeError(`treetime:${channel} returned no bytes`);
    }

    return value;
  }

  function runEvents(id: string, options: TransportEventOptions): Promise<void> {
    const subscriptionId = newSubscriptionId();

    return new Promise<void>((resolve, reject) => {
      const finish = () => {
        ipc.removeListener(RUN_EVENT_CHANNEL, handler);
        options.signal?.removeEventListener("abort", finish);
        resolve();
      };

      const handler = (_event: unknown, eventSubscriptionId: unknown, eventJson: unknown) => {
        if (eventSubscriptionId !== subscriptionId || typeof eventJson !== "string") {
          return;
        }

        try {
          const event = parseRunEvent(JSON.parse(eventJson));

          if (event.type === "log") {
            logToConsole(event.data);
          }

          options.onEvent(event);

          if (event.type === "terminal") {
            finish();
          }
        } catch (error: unknown) {
          ipc.removeListener(RUN_EVENT_CHANNEL, handler);
          reject(error instanceof Error ? error : new Error(String(error)));
        }
      };

      if (options.signal?.aborted === true) {
        resolve();

        return;
      }

      ipc.on(RUN_EVENT_CHANNEL, handler);
      options.signal?.addEventListener("abort", finish);
      ipc.invoke("treetime:runs:subscribe", subscriptionId, id, options.from).catch((error: unknown) => {
        ipc.removeListener(RUN_EVENT_CHANNEL, handler);
        reject(error instanceof Error ? error : new Error(String(error)));
      });
    });
  }

  return {
    version: () => call("version"),
    datasets: () => call("datasets"),
    checkConfig: (request) => call("check-config", JSON.stringify(request)),
    runConfig: (request) => call("run-config", JSON.stringify(request)),
    checkInputs: (request) => call("check-inputs", JSON.stringify(request)),
    listRuns: () => call("runs:list"),
    createRun: (request) => call("runs:create", JSON.stringify(request)),
    getRun: (id) => call("runs:get", id),
    startRun: (id, request) =>
      call("runs:start", id, request.config === undefined ? null : JSON.stringify(request.config)),
    updateRun: (id, request) => call("runs:update", id, JSON.stringify(request)),
    cancelRun: (id) => call("runs:cancel", id),
    deleteRun: async (id) => {
      await call("runs:delete", id);
    },
    restoreRun: (id) => call("runs:restore", id),
    purgeRun: async (id) => {
      await call("runs:purge", id);
    },
    runEvents,
    runFiles: (id) => call("runs:files", id),
    readRunFile: (id, path) => bytes("runs:read-file", id, path),
    runArchive: (id) => bytes("runs:archive", id),
    uploadInput: () => Promise.reject(new LocalInputsError()),
  };
}

function randomSubscriptionId(): string {
  return globalThis.crypto.randomUUID();
}

function decode(value: unknown): unknown {
  if (typeof value !== "string") {
    return value;
  }

  const parsed: unknown = JSON.parse(value);

  return parsed;
}

function logToConsole(log: LogEvent): void {
  switch (log.level) {
    case "error":
      console.error(`[TreeTime] ${log.message}`);
      break;
    case "warn":
      console.warn(`[TreeTime] ${log.message}`);
      break;
    case "info":
    case "debug":
    case "trace":
      console.log(`[TreeTime] [${log.level}] ${log.message}`);
      break;
  }
}
