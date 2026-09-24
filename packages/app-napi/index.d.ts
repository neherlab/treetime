export declare class CommandRunner {
  constructor()
  run(jobId: string, command: string, configJson: string, onEvent: (err: Error | null, eventJson: string) => void): Promise<string>
  cancel(jobId: string): boolean
}

export declare function checkConfigJson(requestJson: string): string

export declare function datasets(): string

export declare function version(): string
