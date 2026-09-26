export declare class Backend {
  constructor(runsDir: string)
  call(requestJson: string): Promise<string>
  subscribe(id: string, from: number, onEvent: (err: Error | null, eventJson: string) => void): Subscription
  fetch(request: PortRequest, onReply: ((arg: PortReply) => void)): PortExchange
  saveRunFile(id: string, path: string, destination: string): Promise<string>
  saveRunArchive(id: string, destination: string): Promise<string>
}

export declare class PortExchange {
  abort(): void
}

export declare class Subscription {
  unsubscribe(): void
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
