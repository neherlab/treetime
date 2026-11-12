import { errorMessage } from "@neherlab/app-contracts";
import { type BlobPortRequest, type PortHeader, type PortReply, blobPortMessageSchema } from "@neherlab/app-ui/port";
import type { MainModule } from "@neherlab/app-wasm";
import { zip } from "remeda";

import { transferList } from "../job/protocol";

export const TRAPPED_MESSAGE = "the back end stopped after an internal error";

export const BODY_SLICE_SIZE = 1 << 20;

export function messagePortEndpoint(port: MessagePort): ServeEndpoint {
  return {
    post: (reply, transfer) => {
      port.postMessage(reply, transfer);
    },
    listen: (listener) => {
      port.addEventListener("message", (event) => {
        listener(event.data);
      });
      port.start();
    },
  };
}

export function responseHeaders({ headerNames, headerValues }: WasmResponseHead): PortHeader[] {
  if (headerNames.length !== headerValues.length) {
    throw new Error(
      `a response head has ${String(headerNames.length)} header names and ${String(headerValues.length)} header values`,
    );
  }

  return zip(headerNames, headerValues).map(([name, value]) => ({ name, value }));
}

export async function* bodyChunks(body: Blob, sliceSize = BODY_SLICE_SIZE): AsyncGenerator<Uint8Array> {
  for (let start = 0; start < body.size; start += sliceSize) {
    yield new Uint8Array(await body.slice(start, start + sliceSize).arrayBuffer());
  }
}

export class ExchangeServer {
  readonly #endpoint: ServeEndpoint;
  readonly #backend: WasmBackend;
  readonly #onTrap: (trap: WebAssembly.RuntimeError) => void;
  readonly #afterHead: () => void;
  readonly #open = new Map<number, WasmExchange>();
  #trapped = false;

  constructor(
    endpoint: ServeEndpoint,
    backend: WasmBackend,
    onTrap: (trap: WebAssembly.RuntimeError) => void,
    afterHead: () => void,
  ) {
    this.#endpoint = endpoint;
    this.#backend = backend;
    this.#onTrap = onTrap;
    this.#afterHead = afterHead;
  }

  listen(): void {
    this.#endpoint.listen((data) => {
      this.#receive(data);
    });
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- port messages arrive untyped; the server parses them with the message schema
  #receive(data: unknown): void {
    const message = blobPortMessageSchema.parse(data);

    if (message.kind === "request") {
      void this.#serve(message.request);
    } else {
      this.#abort(message.seq);
    }
  }

  async #serve(request: BlobPortRequest): Promise<void> {
    const { seq } = request;

    if (this.#trapped) {
      this.#post({ kind: "reset", seq, message: TRAPPED_MESSAGE });

      return;
    }

    let exchange: WasmExchange;

    try {
      // oxlint-disable-next-line custom/require-io-timeout -- the back end runs the request in this worker; the page ends it with an abort message
      exchange = this.#backend.fetch(
        request.method,
        request.url,
        request.headers.map(({ name }) => name),
        request.headers.map(({ value }) => value),
        request.body !== undefined,
      );
    } catch (error) {
      this.#fail(seq, error);

      return;
    }

    this.#open.set(seq, exchange);

    const sending = request.body === undefined ? Promise.resolve(undefined) : this.#send(seq, exchange, request.body);

    try {
      const head = await exchange.head();

      this.#post({ kind: "head", seq, status: head.status, headers: responseHeaders(head) });
      this.#afterHead();

      for (let chunk = await exchange.chunk(); chunk !== undefined; chunk = await exchange.chunk()) {
        this.#post({ kind: "chunk", seq, data: chunk }, transferList(chunk));
      }

      this.#post({ kind: "end", seq });
    } catch (error) {
      const failure = error instanceof WebAssembly.RuntimeError ? error : ((await sending) ?? error);

      if (this.#open.has(seq)) {
        this.#fail(seq, failure);
      }
    }

    this.#open.delete(seq);

    if (!this.#trapped) {
      await sending;
      exchange.free();
    }
  }

  async #send(seq: number, exchange: WasmExchange, body: Blob): Promise<Error | undefined> {
    const error = await sendBody(exchange, body);

    if (error instanceof WebAssembly.RuntimeError && this.#open.has(seq)) {
      this.#open.delete(seq);
      this.#fail(seq, error);
    }

    return error;
  }

  #abort(seq: number): void {
    const exchange = this.#open.get(seq);

    this.#open.delete(seq);
    exchange?.abort();
  }

  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- a rejection carries any value; the error message is read from it
  #fail(seq: number, error: unknown): void {
    if (!(error instanceof WebAssembly.RuntimeError)) {
      this.#post({ kind: "reset", seq, message: errorMessage(error) });

      return;
    }

    this.#post({ kind: "reset", seq, message: TRAPPED_MESSAGE });

    if (!this.#trapped) {
      this.#trapped = true;
      this.#onTrap(error);
    }
  }

  #post(reply: PortReply, transfer: Transferable[] = []): void {
    this.#endpoint.post(reply, transfer);
  }
}

async function sendBody(exchange: WasmExchange, body: Blob): Promise<Error | undefined> {
  try {
    for await (const chunk of bodyChunks(body)) {
      if (!(await exchange.send(chunk))) {
        return undefined;
      }
    }

    exchange.end();

    return undefined;
  } catch (error) {
    const failure = error instanceof Error ? error : new Error(errorMessage(error));

    return failure instanceof WebAssembly.RuntimeError ? failure : abortExchange(exchange, failure);
  }
}

function abortExchange(exchange: WasmExchange, failure: Error): Error {
  try {
    exchange.abort();
  } catch (error) {
    if (error instanceof WebAssembly.RuntimeError) {
      return error;
    }
  }

  return failure;
}

export interface ServeEndpoint {
  post(reply: PortReply, transfer: Transferable[]): void;
  // oxlint-disable-next-line anti-slop/no-unknown-parameters -- port messages arrive untyped; the server parses them
  listen(listener: (data: unknown) => void): void;
}

export interface WasmBackend {
  fetch(...args: Parameters<Backend["fetch"]>): WasmExchange;
}

export interface WasmExchange {
  head(): Promise<WasmResponseHead>;
  chunk: Exchange["chunk"];
  send: Exchange["send"];
  end: Exchange["end"];
  abort: Exchange["abort"];
  free: Exchange["free"];
}

export type WasmResponseHead = Pick<ResponseHead, "status" | "headerNames" | "headerValues">;

type Backend = InstanceType<MainModule["Backend"]>;

type Exchange = ReturnType<Backend["fetch"]>;

type ResponseHead = Awaited<ReturnType<Exchange["head"]>>;
