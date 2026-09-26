import type { PortMessage, PortReply, SaveRunArchiveRequest, SaveRunFileRequest } from "@neherlab/app-napi";

export type SaveRequest =
  | { kind: "save-file"; seq: number; request: SaveRunFileRequest }
  | { kind: "save-archive"; seq: number; request: SaveRunArchiveRequest };

export type ControlRequest = { kind: "port" } | SaveRequest;

export type ControlReply = { kind: "saved"; seq: number } | { kind: "error"; seq: number; error: string };

export interface FetchEndpoint {
  post(message: PortReply): void;
  listen(listener: (message: PortMessage) => void): void;
  onClose(listener: () => void): void;
}
