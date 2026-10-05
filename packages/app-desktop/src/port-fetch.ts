import { ApiError } from "@neherlab/app-contracts/client";
import type { PortHeader, PortMessage, PortReply, PortRequest } from "@neherlab/app-napi";
import type { BackendStopped } from "@neherlab/app-ui/host";

export interface FetchPort {
  postMessage(message: PortMessage): void;
  addEventListener(type: "message", listener: (event: { data: PortReply }) => void): void;
  start(): void;
}

export interface FetchConnection {
  onPort(listener: (port: FetchPort) => void): void;
  onStopped(listener: (stop: BackendStopped) => void): void;
}

export type PortFetch = (input: string | URL | Request, init?: RequestInit) => Promise<Response>;

const SERVER_ERROR_STATUS = 500;

const MISSING_HEAD = "the back end ended the exchange without a response head";

const NULL_BODY_STATUSES = new Set([101, 103, 204, 205, 304]);

export function createPortFetch(connection: FetchConnection): PortFetch {
  const client = new PortClient(connection);

  return async (input, init) => {
    const request = new Request(input, init);
    const body = request.body === null ? undefined : new Uint8Array(await request.arrayBuffer());
    const url = new URL(request.url);

    const portRequest: Omit<PortRequest, "seq"> = {
      method: request.method,
      url: `${url.pathname}${url.search}`,
      headers: Array.from(request.headers, ([name, value]) => ({ name, value })),
    };

    if (body !== undefined) {
      portRequest.body = body;
    }

    return client.exchange(portRequest, request.signal);
  };
}

interface Exchange {
  reply(reply: PortReply): void;
  fail(error: Error): void;
}

class PortClient {
  private port: FetchPort | undefined;
  private readonly waiting = new Map<number, PortRequest>();
  private readonly exchanges = new Map<number, Exchange>();
  private nextSeq = 0;
  private failure: Error | undefined;

  constructor(connection: FetchConnection) {
    connection.onPort((port) => {
      this.connect(port);
    });
    connection.onStopped((stop) => {
      this.stop(stop);
    });
  }

  exchange(request: Omit<PortRequest, "seq">, signal: AbortSignal): Promise<Response> {
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
      this.send({ seq, ...request });
    });
  }

  private send(request: PortRequest): void {
    if (this.port === undefined) {
      this.waiting.set(request.seq, request);
    } else {
      // oxlint-disable-next-line unicorn/require-post-message-target-origin -- a MessagePort takes no target origin
      this.port.postMessage({ kind: "request", request });
    }
  }

  private cancel(seq: number): void {
    if (this.waiting.delete(seq)) {
      return;
    }

    this.port?.postMessage({ kind: "abort", seq });
  }

  private connect(port: FetchPort): void {
    this.port = port;
    this.failure = undefined;
    port.addEventListener("message", (event) => {
      this.exchanges.get(event.data.seq)?.reply(event.data);
    });
    port.start();

    for (const request of this.waiting.values()) {
      port.postMessage({ kind: "request", request });
    }

    this.waiting.clear();
  }

  private stop({ reason, restarts, error: startError }: BackendStopped): void {
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

function abortReason(signal: AbortSignal): Error {
  const reason: unknown = signal.reason;

  return reason instanceof Error ? reason : new DOMException(String(reason), "AbortError");
}
