import type { ErrorResponse } from "@neherlab/app-contracts";
import type { PortMessage, PortReply } from "@neherlab/app-napi";

export type PortScope = "host" | "renderer";

export type ControlRequest = { kind: "port"; scope: PortScope };

export type ControlReply = { kind: "failed"; error: ErrorResponse };

export interface FetchEndpoint {
  post(message: PortReply): void;
  listen(listener: (message: PortMessage) => void): void;
  onClose(listener: () => void): void;
}
