export declare class Backend {
  constructor(runsDir: string)
  call(requestJson: string): Promise<string>
  subscribe(id: string, from: number, onEvent: (err: Error | null, eventJson: string) => void): Subscription
  readFile(id: string, path: string): Buffer
  archive(id: string): Buffer
}

export declare class Subscription {
  unsubscribe(): void
}
