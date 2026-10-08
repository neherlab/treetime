// TypeScript bindings for emscripten-generated code.  Automatically generated at compile time.
interface WasmModule {
}

declare class Backend {
    free(): void;
    [Symbol.dispose](): void;
    fetch(method: string, url: string, headerNames: string[], headerValues: string[], body?: Uint8Array | null): Exchange;
    constructor();
}

declare class Exchange {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    abort(): void;
    chunk(): Promise<Uint8Array | undefined>;
    head(): Promise<ResponseHead>;
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
    setPanicSink(sink: (text: string) => void): void;
    start(): void;
    Backend: typeof Backend;
    Exchange: typeof Exchange;
    ResponseHead: typeof ResponseHead;
}

export { BindgenModule };
export type MainModule = WasmModule & BindgenModule;
export default function MainModuleFactory (options?: unknown): Promise<MainModule>;
