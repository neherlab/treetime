const WASM_FETCH_TIMEOUT_MS = 120_000;

export async function compileWasm(url: string): Promise<WebAssembly.Module> {
  const response = await fetch(url, { signal: AbortSignal.timeout(WASM_FETCH_TIMEOUT_MS) });

  if (!response.ok) {
    throw new Error(`the WebAssembly module answered with HTTP status ${String(response.status)}`);
  }

  return WebAssembly.compileStreaming(response);
}
