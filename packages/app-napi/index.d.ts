export declare class Backend {
  constructor(runsDir: string)
  call(requestJson: string): Promise<string>
  subscribe(id: string, from: number, onEvent: (err: Error | null, eventJson: string) => void): Subscription
  resolveRunFile(id: string, path: string): Promise<string>
  saveRunFile(id: string, path: string, destination: string): Promise<string>
  saveRunArchive(id: string, destination: string): Promise<string>
}

export declare class Subscription {
  unsubscribe(): void
}
