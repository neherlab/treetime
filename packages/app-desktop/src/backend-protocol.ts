import type { ErrorResponse } from "@neherlab/app-contracts";
import type { PortMessage, PortReply, SaveRunArchiveRequest, SaveRunFileRequest } from "@neherlab/app-napi";

export type SaveRequest =
  | { kind: "save-file"; seq: number; request: SaveRunFileRequest }
  | { kind: "save-archive"; seq: number; request: SaveRunArchiveRequest };

export type ControlRequest = { kind: "port" } | SaveRequest;

export type SaveResult = { kind: "saved"; seq: number } | { kind: "error"; seq: number; error: string };

export type ControlReply = SaveResult | { kind: "failed"; error: ErrorResponse };

export interface FetchEndpoint {
  post(message: PortReply): void;
  listen(listener: (message: PortMessage) => void): void;
  onClose(listener: () => void): void;
}
