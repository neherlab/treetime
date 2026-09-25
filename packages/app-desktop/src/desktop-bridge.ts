import {
  createBridge,
  bridgeErrorFromText,
  zPickedFiles,
  type DesktopRequestInput,
  type BridgeTransport,
  type LocalFiles,
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

export const CALL_CHANNEL = "treetime:call";

export type IpcReply = { ok: true; value: unknown } | { ok: false; error: string };

class LocalInputsError extends Error {
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
      return zPickedFiles.parse(unwrap(await ipc.invoke(PICK_FILES_CHANNEL, JSON.stringify(request))));
    },
    pathForFile,
  };
}

function createDesktopTransport(ipc: IpcRendererLike, newSubscriptionId: () => string): BridgeTransport {
  async function call(request: DesktopRequestInput): Promise<unknown> {
    return decode(unwrap(await ipc.invoke(CALL_CHANNEL, JSON.stringify(request))));
  }

  async function bytes(channel: string, ...args: unknown[]): Promise<Uint8Array> {
    const value = unwrap(await ipc.invoke(`treetime:${channel}`, ...args));

    if (!(value instanceof Uint8Array)) {
      throw new TypeError(`treetime:${channel} returned no bytes`);
    }

    return value;
  }

  function runEvents(id: string, options: TransportEventOptions): Promise<void> {
    const subscriptionId = newSubscriptionId();

    return new Promise<void>((resolve, reject) => {
      const stop = () => {
        ipc.removeListener(RUN_EVENT_CHANNEL, handler);
        options.signal?.removeEventListener("abort", finish);
        ipc.send("treetime:runs:unsubscribe", subscriptionId);
      };

      const finish = () => {
        stop();
        resolve();
      };

      const handler = (_event: unknown, eventSubscriptionId: unknown, eventJson: unknown) => {
        if (eventSubscriptionId !== subscriptionId || typeof eventJson !== "string") {
          return;
        }

        try {
          const event: unknown = JSON.parse(eventJson);
          options.onEvent(event);

          if (isTerminal(event)) {
            finish();
          }
        } catch (error: unknown) {
          stop();
          reject(error instanceof Error ? error : new Error(String(error)));
        }
      };

      if (options.signal?.aborted === true) {
        resolve();

        return;
      }

      ipc.on(RUN_EVENT_CHANNEL, handler);
      options.signal?.addEventListener("abort", finish);
      ipc
        .invoke("treetime:runs:subscribe", subscriptionId, id, options.from)
        .then(unwrap)
        .catch((error: unknown) => {
          ipc.removeListener(RUN_EVENT_CHANNEL, handler);
          options.signal?.removeEventListener("abort", finish);
          reject(error instanceof Error ? error : new Error(String(error)));
        });
    });
  }

  return {
    version: () => call({ operation: "version", args: {} }),
    datasets: () => call({ operation: "datasets", args: {} }),
    checkConfig: (request) => call({ operation: "check-config", args: { request } }),
    runConfig: (request) => call({ operation: "run-config", args: { request } }),
    checkInputs: (request) => call({ operation: "check-inputs", args: { request } }),
    listRuns: () => call({ operation: "list-runs", args: {} }),
    createRun: (request) => call({ operation: "create-run", args: { request } }),
    getRun: (id) => call({ operation: "get-run", args: { id } }),
    startRun: (id, request) => call({ operation: "start-run", args: { id, request } }),
    updateRun: (id, request) => call({ operation: "update-run", args: { id, request } }),
    cancelRun: (id) => call({ operation: "cancel-run", args: { id } }),
    deleteRun: async (id) => {
      await call({ operation: "delete-run", args: { id } });
    },
    restoreRun: (id) => call({ operation: "restore-run", args: { id } }),
    purgeRun: async (id) => {
      await call({ operation: "purge-run", args: { id } });
    },
    runEvents,
    runFiles: (id) => call({ operation: "run-files", args: { id } }),
    readRunFile: (id, path) => bytes("runs:read-file", id, path),
    runArchive: (id) => bytes("runs:archive", id),
    uploadInput: () => Promise.reject(new LocalInputsError()),
    runResults: (id) => call({ operation: "run-results", args: { id } }),
    compareRuns: (id, other) => call({ operation: "compare-runs", args: { id, other } }),
    cladeInRuns: (request) => call({ operation: "clade-in-runs", args: { request } }),
  };
}

function randomSubscriptionId(): string {
  return globalThis.crypto.randomUUID();
}

function unwrap(reply: unknown): unknown {
  if (typeof reply !== "object" || reply === null || !("ok" in reply)) {
    throw new TypeError("the main process sent a reply without a result");
  }

  if (reply.ok === true && "value" in reply) {
    return reply.value;
  }

  throw bridgeErrorFromText("error" in reply && typeof reply.error === "string" ? reply.error : "unknown error");
}

function isTerminal(event: unknown): boolean {
  return typeof event === "object" && event !== null && "type" in event && event.type === "terminal";
}

function decode(value: unknown): unknown {
  if (typeof value !== "string") {
    return value;
  }

  const parsed: unknown = JSON.parse(value);

  return parsed;
}
