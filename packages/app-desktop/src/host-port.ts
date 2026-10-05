import type { PortMessage, PortReply } from "@neherlab/app-napi";
import type { BackendStopped } from "@neherlab/app-ui/host";

import type { FetchConnection, FetchPort } from "./port-fetch";

export interface MainPort {
  postMessage(message: PortMessage): void;
  on(event: "message", listener: (event: { data: PortReply }) => void): void;
  start(): void;
}

export function mainFetchPort(port: MainPort): FetchPort {
  return {
    postMessage: (message) => {
      port.postMessage(message);
    },
    addEventListener: (_type, listener) => {
      port.on("message", (event) => {
        listener({ data: event.data });
      });
    },
    start: () => {
      port.start();
    },
  };
}

export class HostConnection implements FetchConnection {
  private readonly portListeners: Array<(port: FetchPort) => void> = [];
  private readonly stopListeners: Array<(stop: BackendStopped) => void> = [];
  private current: FetchPort | undefined;

  connect(port: FetchPort): void {
    this.current = port;
    this.portListeners.forEach((listener) => {
      listener(port);
    });
  }

  stop(stop: BackendStopped): void {
    this.current = undefined;
    this.stopListeners.forEach((listener) => {
      listener(stop);
    });
  }

  onPort(listener: (port: FetchPort) => void): void {
    this.portListeners.push(listener);

    if (this.current !== undefined) {
      listener(this.current);
    }
  }

  onStopped(listener: (stop: BackendStopped) => void): void {
    this.stopListeners.push(listener);
  }
}
