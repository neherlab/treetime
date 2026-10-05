export declare class Backend {
  constructor()
  fetch(request: PortRequest, scope: 'host' | 'renderer', onReply: ((arg: PortReply) => void)): PortExchange
  rejectRequest(seq: number, message: string): Array<PortReply>
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
  | { kind: 'reset'; seq: number; message: string }

export interface PortRequest {
  seq: number
  method: string
  url: string
  headers: Array<PortHeader>
  body?: Uint8Array
}
