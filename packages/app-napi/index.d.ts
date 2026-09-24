export declare function ancestral(argsJson: string, onEvent: (err: Error | null, eventJson: string) => void): Promise<string>

export declare function ancestralSync(argsJson: string): string

export declare function cancel(): void

export declare function clock(argsJson: string, onEvent: (err: Error | null, eventJson: string) => void): Promise<string>

export declare function datasets(): string

export declare function mugration(argsJson: string, onEvent: (err: Error | null, eventJson: string) => void): Promise<string>

export declare function optimize(argsJson: string, onEvent: (err: Error | null, eventJson: string) => void): Promise<string>

export declare function prune(argsJson: string, onEvent: (err: Error | null, eventJson: string) => void): Promise<string>

export declare function timetree(argsJson: string, onEvent: (err: Error | null, eventJson: string) => void): Promise<string>

export declare function version(): string
