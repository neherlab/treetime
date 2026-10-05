export declare class Backend {
  constructor()
  fetch(request: PortRequest, onReply: ((arg: PortReply) => void)): PortExchange
  saveRunFile(request: SaveRunFileRequest): Promise<void>
  saveRunArchive(request: SaveRunArchiveRequest): Promise<void>
}

export declare class PortExchange {
  abort(): void
}

export declare function appStartup(): AppStartup

export interface AppStartup {
  profileDir: string
  logsDir: string
  theme: string
}

export interface PortError {
  code: string
  message: string
  causes: Array<string>
}

export interface PortHeader {
  name: string
  value: string
}

export type PortMessage =
  | { kind: 'request'; request: PortRequest }
  | { kind: 'abort'; seq: number }

export type PortReply =
  | { kind: 'head'; seq: number; status: number; headers: Array<PortHeader> }
  | { kind: 'chunk'; seq: number; data: Uint8Array }
  | { kind: 'end'; seq: number }
  | { kind: 'error'; seq: number; error: PortError }

export interface PortRequest {
  seq: number
  method: string
  url: string
  headers: Array<PortHeader>
  body?: string
}

export interface SaveRunArchiveRequest {
  id: string
  destination: string
}

export interface SaveRunFileRequest {
  id: string
  path: string
  destination: string
}
