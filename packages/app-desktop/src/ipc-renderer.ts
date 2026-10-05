import {
  BACKEND_PORT_CHANNEL,
  BACKEND_PORT_REQUEST_CHANNEL,
  hostChannelName,
  hostChannelSchemas,
  hostEventSchema,
  type BackendStopped,
  type Host,
  type HostChannel,
  type HostReply,
  type HostRequest,
} from "@neherlab/app-ui/host";

import type { FetchConnection, FetchPort } from "./port-fetch";

export interface IpcRendererLike<Port> {
  // oxlint-disable-next-line anti-slop/no-unknown-returns -- Electron IPC replies are untyped; invoke() parses them with the channel schema
  invoke(channel: string, ...args: unknown[]): Promise<unknown>;
  send(channel: string, ...args: unknown[]): void;
  on(channel: string, listener: (event: { ports: Port[] }, ...args: unknown[]) => void): void;
}

export interface FileUtils {
  getPathForFile(file: File): string;
}

export interface PortTarget<Port> {
  postMessage(message: { channel: string }, targetOrigin: string, transfer: Port[]): void;
}

interface WindowMessage {
  source: unknown;
  data: unknown;
  ports: readonly [FetchPort?];
}

export interface WindowLike {
  addEventListener(type: "message", listener: (event: WindowMessage) => void): void;
}

export async function invoke<Port, C extends HostChannel>(
  ipc: IpcRendererLike<Port>,
  channel: C,
  request: HostRequest<C>,
): Promise<HostReply<C>> {
  const schemas = hostChannelSchemas(channel);

  return schemas.reply.parse(await ipc.invoke(hostChannelName(channel), schemas.request.parse(request)));
}

export function subscribeBackendStopped<Port>(
  ipc: IpcRendererLike<Port>,
  listener: (stop: BackendStopped) => void,
): void {
  const schema = hostEventSchema("backend-stopped");

  ipc.on(hostChannelName("backend-stopped"), (_event, payload) => {
    listener(schema.parse(payload));
  });
}

export function relayBackendPort<Port>(ipc: IpcRendererLike<Port>, target: PortTarget<Port>): void {
  ipc.on(BACKEND_PORT_CHANNEL, (event) => {
    target.postMessage({ channel: BACKEND_PORT_CHANNEL }, "*", event.ports);
  });
}

export function createPreloadHost<Port>(ipc: IpcRendererLike<Port>, files: FileUtils): Host {
  return {
    pickFiles: async (request) => invoke(ipc, "pick-files", request),
    pickFolder: async (request) => invoke(ipc, "pick-folder", request),
    saveRun: async (request) => invoke(ipc, "save-run", request),
    restartBackend: async () => invoke(ipc, "restart-backend", undefined),
    setNativeTheme: async (theme) => invoke(ipc, "set-native-theme", theme),
    pathForFile: (file) => files.getPathForFile(file),
    connectBackend: () => {
      ipc.send(BACKEND_PORT_REQUEST_CHANNEL);
    },
    onBackendStopped: (listener) => {
      subscribeBackendStopped(ipc, listener);
    },
  };
}

export function windowFetchConnection(
  target: WindowLike,
  host: Pick<Host, "connectBackend" | "onBackendStopped">,
): FetchConnection {
  const portListeners: Array<(port: FetchPort) => void> = [];
  let current: FetchPort | undefined;
  let requested = false;

  host.onBackendStopped(() => {
    current = undefined;
  });

  target.addEventListener("message", (event) => {
    const [port] = event.ports;

    if (event.source !== target || !isPortMessage(event.data) || port === undefined) {
      return;
    }

    current = port;
    portListeners.forEach((listener) => {
      listener(port);
    });
  });

  return {
    onPort(listener) {
      portListeners.push(listener);

      if (current !== undefined) {
        listener(current);
      } else if (!requested) {
        requested = true;
        host.connectBackend();
      }
    },
    onStopped(listener) {
      host.onBackendStopped(listener);
    },
  };
}

// oxlint-disable-next-line anti-slop/no-unknown-parameters -- window messages carry untyped data; this checks the shape of the port message
function isPortMessage(data: unknown): boolean {
  // oxlint-disable-next-line anti-slop/no-runtime-typeof -- the port message is a plain marker object posted by the preload script
  return typeof data === "object" && data !== null && "channel" in data && data.channel === BACKEND_PORT_CHANNEL;
}
