import { ApiError } from "@neherlab/app-contracts/client";

import type { BackendStopped } from "./host";
import { type BlobPortRequest, type PortHeader, type PortReply, type PortRequest, portReplySchema } from "./port";

export interface FetchPort<R extends BlobPortRequest = PortRequest> {
  postMessage(message: PortMessageOf<R>): void;
  addEventListener(type: "message", listener: (event: { data: unknown }) => void): void;
  start(): void;
}

export interface FetchConnection<R extends BlobPortRequest = PortRequest> {
  onPort(listener: (port: FetchPort<R>) => void): void;
  onStopped(listener: (stop: BackendStopped) => void): void;
}

export type PortMessageOf<R extends BlobPortRequest> = { kind: "request"; request: R } | { kind: "abort"; seq: number };

export type PortFetch = (input: string | URL | Request, init?: RequestInit) => Promise<Response>;

export type StopListener = (stop: BackendStopped) => void;

const SERVER_ERROR_STATUS = 500;

const MISSING_HEAD = "the back end ended the exchange without a response head";

const MALFORMED_REPLY = "the back end sent a malformed message";

const NULL_BODY_STATUSES = new Set([101, 103, 204, 205, 304]);

export function createPortFetch(connection: FetchConnection, onStop: StopListener = ignoreStop): PortFetch {
  const client = new PortClient(connection, onStop);

  return async (input, init) => {
    const request = new Request(input, init);
    const head = requestHead(request);

    const prepared: Omit<PortRequest, "seq"> =
      request.body === null ? head : { ...head, body: new Uint8Array(await request.arrayBuffer()) };

    return client.exchange((seq) => ({ seq, ...prepared }), request.signal);
  };
}

export function createBlobPortFetch(
  connection: FetchConnection<BlobPortRequest>,
  onStop: StopListener = ignoreStop,
): PortFetch {
  const client = new PortClient(connection, onStop);

  return async (input, init) => {
    const request = new Request(input, init);
    const head = requestHead(request);

    const prepared: Omit<BlobPortRequest, "seq"> =
      request.body === null ? head : { ...head, body: await request.blob() };

    return client.exchange((seq) => ({ seq, ...prepared }), request.signal);
  };
}

function requestHead(request: Request): Omit<PortRequest, "seq" | "body"> {
  const url = new URL(request.url);

  return {
    method: request.method,
    url: `${url.pathname}${url.search}`,
    headers: Array.from(request.headers, ([name, value]) => ({ name, value })),
  };
}

interface Exchange {
  reply(reply: PortReply): void;
  fail(error: Error): void;
}

class PortClient<R extends BlobPortRequest> {
  private port: FetchPort<R> | undefined;
  private readonly waiting = new Map<number, R>();
  private readonly exchanges = new Map<number, Exchange>();
  private nextSeq = 0;
  private failure: Error | undefined;
  private readonly onStop: StopListener;

  constructor(connection: FetchConnection<R>, onStop: StopListener) {
    this.onStop = onStop;
    connection.onPort((port) => {
      this.connect(port);
    });
    connection.onStopped((stop) => {
      this.stop(stop);
    });
  }

  exchange(build: (seq: number) => R, signal: AbortSignal): Promise<Response> {
    if (signal.aborted) {
      return Promise.reject(abortReason(signal));
    }

    if (this.failure !== undefined) {
      return Promise.reject(this.failure);
    }

    const seq = this.nextSeq;
    this.nextSeq += 1;

    return new Promise((resolve, reject) => {
      let responded = false;
      let body: ReadableStreamDefaultController<Uint8Array> | undefined;

      const close = () => {
        this.exchanges.delete(seq);
        signal.removeEventListener("abort", abort);
      };

      const fail = (error: Error) => {
        close();

        if (responded) {
          body?.error(error);
        } else {
          reject(error);
        }
      };

      const abort = () => {
        this.cancel(seq);
        fail(abortReason(signal));
      };

      this.exchanges.set(seq, {
        reply: (reply) => {
          switch (reply.kind) {
            case "head":
              responded = true;
              resolve(
                response(
                  reply.status,
                  reply.headers,
                  (controller) => {
                    body = controller;
                  },
                  () => {
                    this.cancel(seq);
                    close();
                  },
                ),
              );
              break;
            case "chunk":
              body?.enqueue(reply.data);
              break;
            case "end":
              if (responded) {
                close();
                body?.close();
              } else {
                fail(new TypeError(MISSING_HEAD));
              }

              break;
            case "reset":
              fail(new TypeError(reply.message));
              break;
          }
        },
        fail,
      });
      signal.addEventListener("abort", abort, { once: true });
      this.send(build(seq));
    });
  }

  private send(request: R): void {
    if (this.port === undefined) {
      this.waiting.set(request.seq, request);
    } else {
      this.post(this.port, request);
    }
  }

  private cancel(seq: number): void {
    if (this.waiting.delete(seq)) {
      return;
    }

    this.port?.postMessage({ kind: "abort", seq });
  }

  private post(port: FetchPort<R>, request: R): void {
    port.postMessage({ kind: "request", request });
  }

  private connect(port: FetchPort<R>): void {
    this.port = port;
    this.failure = undefined;
    port.addEventListener("message", (event) => {
      const reply = portReplySchema.safeParse(event.data);

      if (reply.success) {
        this.exchanges.get(reply.data.seq)?.reply(reply.data);
      } else {
        this.fault(new TypeError(MALFORMED_REPLY, { cause: reply.error }));
      }
    });
    port.start();

    for (const request of this.waiting.values()) {
      this.post(port, request);
    }

    this.waiting.clear();
  }

  private fault(error: Error): void {
    this.port = undefined;
    this.failure = error;
    this.waiting.clear();

    for (const exchange of this.exchanges.values()) {
      exchange.fail(error);
    }
  }

  private stop(stop: BackendStopped): void {
    const { reason, restarts, error: startError } = stop;
    this.port = undefined;
    const error = startError === undefined ? new TypeError(reason) : new ApiError(SERVER_ERROR_STATUS, startError);

    if (!restarts) {
      this.failure = error;
      this.waiting.clear();
    }

    for (const [seq, exchange] of this.exchanges) {
      if (!this.waiting.has(seq)) {
        exchange.fail(error);
      }
    }

    this.onStop(stop);
  }
}

function response(
  status: number,
  headers: readonly PortHeader[],
  start: (controller: ReadableStreamDefaultController<Uint8Array>) => void,
  cancel: () => void,
): Response {
  const init = { status, headers: headers.map(({ name, value }): [string, string] => [name, value]) };

  if (NULL_BODY_STATUSES.has(status)) {
    return new Response(null, init);
  }

  return new Response(new ReadableStream<Uint8Array>({ start, cancel }), init);
}

function ignoreStop(): void {}

function abortReason(signal: AbortSignal): Error {
  const reason: unknown = signal.reason;

  return reason instanceof Error ? reason : new DOMException(String(reason), "AbortError");
}
