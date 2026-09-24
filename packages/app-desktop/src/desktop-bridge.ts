import {
  CancelledError,
  createBridge,
  parseJobEvent,
  type AppCommand,
  type BridgeTransport,
  type LogEvent,
  type TransportCommandOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";

export interface IpcRendererLike {
  invoke(channel: string, ...args: unknown[]): Promise<unknown>;
  on(channel: string, listener: (event: unknown, ...args: unknown[]) => void): void;
  removeListener(channel: string, listener: (event: unknown, ...args: unknown[]) => void): void;
  send(channel: string, ...args: unknown[]): void;
}

export const JOB_EVENT_CHANNEL = "treetime:job-event";

export function createDesktopBridge(ipc: IpcRendererLike, newJobId: () => string = randomJobId): TreeTimeBridge {
  return createBridge(createDesktopTransport(ipc, newJobId));
}

function createDesktopTransport(ipc: IpcRendererLike, newJobId: () => string): BridgeTransport {
  async function query(endpoint: string): Promise<unknown> {
    return decode(await ipc.invoke(`treetime:${endpoint}`));
  }

  async function request(endpoint: string, body: unknown): Promise<unknown> {
    return decode(await ipc.invoke(`treetime:${endpoint}`, JSON.stringify(body)));
  }

  async function command(command: AppCommand, config: unknown, options: TransportCommandOptions): Promise<unknown> {
    const jobId = newJobId();

    const eventHandler = (_event: unknown, eventJobId: unknown, eventJson: unknown) => {
      if (eventJobId !== jobId || typeof eventJson !== "string") {
        return;
      }

      const event = parseJobEvent(JSON.parse(eventJson));

      if (event.type === "log") {
        logToConsole(event.data);
      }

      options.onEvent(event);
    };

    const abortHandler = () => {
      ipc.send("treetime:cancel", jobId);
    };

    if (options.signal?.aborted === true) {
      throw new CancelledError();
    }

    ipc.on(JOB_EVENT_CHANNEL, eventHandler);
    options.signal?.addEventListener("abort", abortHandler);

    try {
      return decode(await ipc.invoke("treetime:run", jobId, command, JSON.stringify(config)));
    } finally {
      ipc.removeListener(JOB_EVENT_CHANNEL, eventHandler);
      options.signal?.removeEventListener("abort", abortHandler);
    }
  }

  return { query, request, command };
}

function randomJobId(): string {
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
