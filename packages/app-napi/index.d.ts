export declare class RunService {
  constructor(runsDir: string)
  list(): string
  get(id: string): string
  create(requestJson: string): string
  start(id: string, configJson: string | null): Promise<string>
  cancel(id: string): boolean
  subscribe(id: string, from: number, onEvent: (err: Error | null, eventJson: string) => void): Subscription
  update(id: string, requestJson: string): string
  delete(id: string): void
  restore(id: string): string
  purge(id: string): void
  files(id: string): string
  readFile(id: string, path: string): Buffer
  archive(id: string): Buffer
  results(id: string): Promise<string>
  compare(id: string, other: string): Promise<string>
  cladeInRuns(requestJson: string): Promise<string>
}

export declare class Subscription {
  unsubscribe(): void
}

export declare function checkConfigJson(requestJson: string): string

export declare function checkInputsJson(requestJson: string): Promise<string>

export declare function datasets(): string

export declare function runConfigJson(requestJson: string): string

export declare function version(): string
