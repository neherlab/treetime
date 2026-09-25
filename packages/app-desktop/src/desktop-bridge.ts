import {
  BridgeError,
  bridgeErrorFromText,
  createBridge,
  zPickedFiles,
  type BridgeTransport,
  type DesktopRequestInput,
  type LocalFiles,
  type PickFilesRequest,
  type TransportEventOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";
import * as z from "zod";

import {
  portEndpoint,
  zBackendReply,
  type BackendReply,
  type BackendRequest,
  type ClientEndpoint,
  type PortLike,
} from "./backend-protocol";
import { BACKEND_PORT_CHANNEL, BACKEND_STOPPED_CHANNEL } from "./channels";

export const zShellMessage = z.discriminatedUnion("channel", [
  z.object({ channel: z.literal(BACKEND_PORT_CHANNEL) }),
  z.object({ channel: z.literal(BACKEND_STOPPED_CHANNEL), reason: z.string() }),
]);

export type ShellMessage = z.infer<typeof zShellMessage>;

export interface BackendConnection {
  onEndpoint(listener: (endpoint: ClientEndpoint) => void): void;
  onStopped(listener: (reason: string) => void): void;
}

export interface WindowMessage {
  source: unknown;
  data: unknown;
  ports: readonly PortLike[];
}

export interface WindowLike {
  addEventListener(type: "message", listener: (event: WindowMessage) => void): void;
}

export interface DesktopShell {
  connectBackend(): void;
  pickFiles(request: PickFilesRequest): Promise<unknown>;
  pathForFile(file: File): string;
}

class LocalInputsError extends Error {
  constructor() {
    super("the desktop application reads inputs from local file paths; name the files in the run configuration");
    this.name = "LocalInputsError";
  }
}

export function createDesktopBridge(connection: BackendConnection): TreeTimeBridge {
  return createBridge(createPortTransport(new BackendClient(connection)));
}

export function createLocalFiles(shell: DesktopShell): LocalFiles {
  return {
    async pickFiles(request) {
      return zPickedFiles.parse(await shell.pickFiles(request));
    },
    pathForFile: (file) => shell.pathForFile(file),
  };
}

export function windowBackendConnection(target: WindowLike, shell: DesktopShell): BackendConnection {
  const endpointListeners: Array<(endpoint: ClientEndpoint) => void> = [];
  const stoppedListeners: Array<(reason: string) => void> = [];

  target.addEventListener("message", (event) => {
    const message = zShellMessage.safeParse(event.data);
    const [port] = event.ports;

    if (event.source !== target || !message.success) {
      return;
    }

    if (message.data.channel === BACKEND_STOPPED_CHANNEL) {
      const { reason } = message.data;
      stoppedListeners.forEach((listener) => {
        listener(reason);
      });
    } else if (port !== undefined) {
      const endpoint = portEndpoint<BackendReply, BackendRequest>(port, zBackendReply);
      endpointListeners.forEach((listener) => {
        listener(endpoint);
      });
    }
  });

  return {
    onEndpoint(listener) {
      endpointListeners.push(listener);
      shell.connectBackend();
    },
    onStopped(listener) {
      stoppedListeners.push(listener);
    },
  };
}

interface Handler {
  reply(reply: BackendReply): void;
  stopped(error: BridgeError): void;
  resume?: () => BackendRequest;
}

class BackendClient {
  private endpoint: ClientEndpoint | undefined;
  private readonly waiting: BackendRequest[] = [];
  private readonly handlers = new Map<number, Handler>();
  private nextSeq = 0;

  constructor(connection: BackendConnection) {
    connection.onEndpoint((endpoint) => {
      this.connect(endpoint);
    });
    connection.onStopped((reason) => {
      this.stop(reason);
    });
  }

  open(handler: Handler): number {
    const seq = this.nextSeq;
    this.nextSeq += 1;
    this.handlers.set(seq, handler);

    return seq;
  }

  close(seq: number): void {
    this.handlers.delete(seq);
  }

  send(request: BackendRequest): void {
    if (this.endpoint === undefined) {
      this.waiting.push(request);
    } else {
      this.endpoint.post(request);
    }
  }

  private connect(endpoint: ClientEndpoint): void {
    this.endpoint = endpoint;
    endpoint.listen((reply) => {
      this.handlers.get(reply.seq)?.reply(reply);
    });

    for (const handler of this.handlers.values()) {
      if (handler.resume !== undefined) {
        endpoint.post(handler.resume());
      }
    }

    for (const request of this.waiting.splice(0)) {
      endpoint.post(request);
    }
  }

  private stop(reason: string): void {
    this.endpoint = undefined;
    const error = new BridgeError({ code: "internal_error", message: reason, causes: [] });

    for (const [seq, handler] of this.handlers) {
      if (handler.resume === undefined) {
        this.handlers.delete(seq);
        handler.stopped(error);
      }
    }
  }
}

function createPortTransport(client: BackendClient): BridgeTransport {
  function call(request: DesktopRequestInput): Promise<unknown> {
    return new Promise((resolve, reject) => {
      const seq = client.open({
        reply(reply) {
          client.close(seq);

          if (reply.kind === "result") {
            resolve(JSON.parse(reply.json));
          } else if (reply.kind === "error") {
            reject(bridgeErrorFromText(reply.error));
          } else {
            reject(new TypeError(`the back end answered a call with a ${reply.kind} message`));
          }
        },
        stopped: reject,
      });

      client.send({ kind: "call", seq, request: JSON.stringify(request) });
    });
  }

  function bytes(request: (seq: number) => BackendRequest): Promise<Uint8Array> {
    return new Promise((resolve, reject) => {
      const chunks: Uint8Array[] = [];

      const seq = client.open({
        reply(reply) {
          if (reply.kind === "chunk") {
            chunks.push(new Uint8Array(reply.bytes));

            return;
          }

          client.close(seq);

          if (reply.kind === "end") {
            resolve(concat(chunks));
          } else if (reply.kind === "error") {
            reject(bridgeErrorFromText(reply.error));
          } else {
            reject(new TypeError(`the back end answered a file read with a ${reply.kind} message`));
          }
        },
        stopped: reject,
      });

      client.send(request(seq));
    });
  }

  function runEvents(id: string, options: TransportEventOptions): Promise<void> {
    return new Promise((resolve, reject) => {
      if (options.signal?.aborted === true) {
        resolve();

        return;
      }

      let from = options.from;

      const end = () => {
        client.close(seq);
        client.send({ kind: "unsubscribe", seq });
        options.signal?.removeEventListener("abort", finish);
      };

      const finish = () => {
        end();
        resolve();
      };

      const seq = client.open({
        reply(reply) {
          if (reply.kind === "error") {
            end();
            reject(bridgeErrorFromText(reply.error));

            return;
          }

          if (reply.kind !== "event") {
            return;
          }

          try {
            const event = options.onEvent(JSON.parse(reply.json));
            from = event.seq + 1;

            if (event.type === "terminal") {
              finish();
            }
          } catch (error: unknown) {
            end();
            reject(error instanceof Error ? error : new Error(String(error)));
          }
        },
        stopped: reject,
        resume: () => ({ kind: "subscribe", seq, id, from }),
      });

      options.signal?.addEventListener("abort", finish);
      client.send({ kind: "subscribe", seq, id, from });
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
    readRunFile: (id, path) => bytes((seq) => ({ kind: "read-file", seq, id, path })),
    runArchive: (id) => bytes((seq) => ({ kind: "archive", seq, id })),
    uploadInput: () => Promise.reject(new LocalInputsError()),
    runResults: (id) => call({ operation: "run-results", args: { id } }),
    compareRuns: (id, other) => call({ operation: "compare-runs", args: { id, other } }),
    cladeInRuns: (request) => call({ operation: "clade-in-runs", args: { request } }),
  };
}

function concat(chunks: Uint8Array[]): Uint8Array {
  const result = new Uint8Array(chunks.reduce((size, chunk) => size + chunk.byteLength, 0));
  let offset = 0;

  for (const chunk of chunks) {
    result.set(chunk, offset);
    offset += chunk.byteLength;
  }

  return result;
}
