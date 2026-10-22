// TypeScript bindings for emscripten-generated code.  Automatically generated at compile time.
interface WasmModule {
}

declare class Backend {
    free(): void;
    [Symbol.dispose](): void;
    appendRunEvent(runId: string, eventJson: string): void;
    appendRunOutput(runId: string, relativePath: string, chunk: Uint8Array): void;
    beginRunOutput(runId: string, relativePath: string): void;
    completeRun(runId: string, outcomeJson: string): boolean;
    failRun(runId: string, lines: string[]): boolean;
    fetch(method: string, url: string, headerNames: string[], headerValues: string[], body?: Uint8Array | null): Exchange;
    constructor();
    readRunInput(runId: string, path: string): ChunkReader;
    takeJobOrders(): JobOrder[];
}

declare class ChunkReader {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    chunk(): Uint8Array | undefined;
}

declare class Exchange {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    abort(): void;
    chunk(): Promise<Uint8Array | undefined>;
    head(): Promise<ResponseHead>;
}

declare class Job {
    free(): void;
    [Symbol.dispose](): void;
    appendInput(path: string, chunk: Uint8Array): void;
    beginInput(path: string): void;
    constructor(jobJson: string, jobs: number);
    outputFiles(): string[];
    readOutput(path: string): ChunkReader;
    run(onEvent: (event: string) => void): Promise<string>;
}

declare class JobOrder {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    inputPaths: string[];
    get job(): string | undefined;
    set job(value: string | null | undefined);
    kind: string;
    runId: string;
}

declare class ResponseHead {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    headerNames: string[];
    headerValues: string[];
    status: number;
}

interface BindgenModule {
    /* tslint:disable */
    /* eslint-disable */
    start(): void;
    Backend: typeof Backend;
    ChunkReader: typeof ChunkReader;
    Exchange: typeof Exchange;
    Job: typeof Job;
    JobOrder: typeof JobOrder;
    ResponseHead: typeof ResponseHead;
}

export { BindgenModule };
export type MainModule = WasmModule & BindgenModule;
export default function MainModuleFactory (options?: unknown): Promise<MainModule>;
