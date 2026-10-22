import { errorMessage } from "@neherlab/app-contracts";
import { type PortHeader, type PortReply, type PortRequest, portMessageSchema } from "@neherlab/app-ui/port";
import type { MainModule } from "@neherlab/app-wasm";
import { zip } from "remeda";

import { transferList } from "../job/protocol";

export const TRAPPED_MESSAGE = "the back end stopped after an internal error";

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
    const message = portMessageSchema.parse(data);

    if (message.kind === "request") {
      void this.#serve(message.request);
    } else {
      this.#abort(message.seq);
    }
  }

  async #serve(request: PortRequest): Promise<void> {
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
        request.body,
      );
    } catch (error) {
      this.#fail(seq, error);

      return;
    }

    this.#open.set(seq, exchange);

    try {
      const head = await exchange.head();

      this.#post({ kind: "head", seq, status: head.status, headers: responseHeaders(head) });
      this.#afterHead();

      for (let chunk = await exchange.chunk(); chunk !== undefined; chunk = await exchange.chunk()) {
        this.#post({ kind: "chunk", seq, data: chunk }, transferList(chunk));
      }

      this.#post({ kind: "end", seq });
    } catch (error) {
      if (this.#open.has(seq)) {
        this.#fail(seq, error);
      }
    } finally {
      this.#open.delete(seq);
      exchange.free();
    }
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
  abort: Exchange["abort"];
  free: Exchange["free"];
}

export type WasmResponseHead = Pick<ResponseHead, "status" | "headerNames" | "headerValues">;

type Backend = InstanceType<MainModule["Backend"]>;

type Exchange = ReturnType<Backend["fetch"]>;

type ResponseHead = Awaited<ReturnType<Exchange["head"]>>;
