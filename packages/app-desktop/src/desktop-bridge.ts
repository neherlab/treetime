import {
  CancelledError,
  createBridge,
  parseLogEvent,
  parseProgressEvent,
  type BridgeTransport,
  type CommandOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";

export interface IpcRendererLike {
  invoke(channel: string, ...args: unknown[]): Promise<unknown>;
  on(channel: string, listener: (event: unknown, ...args: unknown[]) => void): void;
  removeListener(channel: string, listener: (event: unknown, ...args: unknown[]) => void): void;
  send(channel: string, ...args: unknown[]): void;
}

export function createDesktopBridge(ipc: IpcRendererLike): TreeTimeBridge {
  return createBridge(createDesktopTransport(ipc));
}

export function createDesktopTransport(ipc: IpcRendererLike): BridgeTransport {
  async function query(endpoint: string): Promise<unknown> {
    return decode(await ipc.invoke(`treetime:${endpoint}`));
  }

  async function command(endpoint: string, args: unknown, options?: CommandOptions): Promise<unknown> {
    const progressHandler = (_event: unknown, data: unknown) => {
      options?.onProgress?.(parseProgressEvent(data));
    };

    const logHandler = (_event: unknown, data: unknown) => {
      logToConsole(parseLogEvent(data));
    };

    const abortHandler = () => {
      ipc.send("treetime:cancel");
    };

    ipc.on("treetime:progress", progressHandler);
    ipc.on("treetime:log", logHandler);
    options?.signal?.addEventListener("abort", abortHandler);

    try {
      return decode(await ipc.invoke(`treetime:${endpoint}`, JSON.stringify(args)));
    } catch (err: unknown) {
      if (err instanceof Error && err.message.includes("cancelled")) {
        throw new CancelledError();
      }

      throw err;
    } finally {
      ipc.removeListener("treetime:progress", progressHandler);
      ipc.removeListener("treetime:log", logHandler);
      options?.signal?.removeEventListener("abort", abortHandler);
    }
  }

  return { query, command };
}

function decode(value: unknown): unknown {
  if (typeof value !== "string") {
    return value;
  }

  const parsed: unknown = JSON.parse(value);

  return parsed;
}

function logToConsole(log: { level: string; message: string }): void {
  switch (log.level) {
    case "Error":
      console.error(`[TreeTime] ${log.message}`);
      break;
    case "Warn":
      console.warn(`[TreeTime] ${log.message}`);
      break;
    default:
      console.log(`[TreeTime] [${log.level}] ${log.message}`);
      break;
  }
}
