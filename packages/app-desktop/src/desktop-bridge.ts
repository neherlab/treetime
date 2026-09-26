import {
  BridgeError,
  bridgeErrorFromText,
  createBridge,
  zPickedFiles,
  type BridgeTransport,
  type LocalFiles,
  type OperationRequestInput,
  type PickFilesRequest,
  type TransportEventOptions,
  type TreeTimeBridge,
} from "@neherlab/app-contracts";

import {
  portEndpoint,
  zBackendReply,
  type BackendReply,
  type BackendRequest,
  type ClientEndpoint,
  type PortLike,
} from "./backend-protocol";
import { BACKEND_STOPPED_CHANNEL } from "./channels";
import {
  zSaveReply,
  zShellMessage,
  type SaveReply,
  type SaveRunArchiveRequest,
  type SaveRunFileRequest,
} from "./shell-protocol";

export interface BackendConnection {
  onEndpoint(listener: (endpoint: ClientEndpoint) => void): void;
  onStopped(listener: (reason: string, restarts: boolean) => void): void;
}

interface WindowMessage {
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
  saveRunFile(request: SaveRunFileRequest): Promise<unknown>;
  saveRunArchive(request: SaveRunArchiveRequest): Promise<unknown>;
}

class LocalInputsError extends Error {
  constructor() {
    super("the desktop application reads inputs from local file paths; name the files in the run configuration");
    this.name = "LocalInputsError";
  }
}

export function createDesktopBridge(connection: BackendConnection, shell: DesktopShell): TreeTimeBridge {
  return createBridge(createPortTransport(new BackendClient(connection), shell));
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
  const stoppedListeners: Array<(reason: string, restarts: boolean) => void> = [];

  target.addEventListener("message", (event) => {
    const message = zShellMessage.safeParse(event.data);
    const [port] = event.ports;

    if (event.source !== target || !message.success) {
      return;
    }

    if (message.data.channel === BACKEND_STOPPED_CHANNEL) {
      const { reason, restarts } = message.data;
      stoppedListeners.forEach((listener) => {
        listener(reason, restarts);
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
  private failure: BridgeError | undefined;

  constructor(connection: BackendConnection) {
    connection.onEndpoint((endpoint) => {
      this.connect(endpoint);
    });
    connection.onStopped((reason, restarts) => {
      this.stop(reason, restarts);
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
    if (this.failure !== undefined) {
      const handler = this.handlers.get(request.seq);
      this.handlers.delete(request.seq);
      handler?.stopped(this.failure);
    } else if (this.endpoint === undefined) {
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

  private stop(reason: string, restarts: boolean): void {
    this.endpoint = undefined;
    const error = new BridgeError({ code: "internal_error", message: reason, causes: [] });

    if (!restarts) {
      this.failure = error;
      this.waiting.splice(0);
    }

    for (const [seq, handler] of this.handlers) {
      if (!restarts || handler.resume === undefined) {
        this.handlers.delete(seq);
        handler.stopped(error);
      }
    }
  }
}

function createPortTransport(client: BackendClient, shell: DesktopShell): BridgeTransport {
  function call(request: OperationRequestInput): Promise<unknown> {
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
    call,
    runEvents,
    saveRunFile: async (id, path, name) => saved(zSaveReply.parse(await shell.saveRunFile({ id, path, name }))),
    saveRunArchive: async (id, name) => saved(zSaveReply.parse(await shell.saveRunArchive({ id, name }))),
    uploadInput: () => Promise.reject(new LocalInputsError()),
  };
}

function saved(reply: SaveReply): boolean {
  if ("error" in reply) {
    throw bridgeErrorFromText(reply.error);
  }

  return reply.saved;
}
