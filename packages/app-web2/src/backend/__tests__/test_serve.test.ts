import type { PortReply } from "@neherlab/app-ui/port";
import pWaitFor from "p-wait-for";

import {
  ExchangeServer,
  responseHeaders,
  TRAPPED_MESSAGE,
  type ServeEndpoint,
  type WasmBackend,
  type WasmExchange,
  type WasmResponseHead,
} from "../serve";

const REQUEST = { kind: "request", request: { seq: 7, method: "GET", url: "/api/version", headers: [] } };

describe("serve exchanges", () => {
  test("a request goes to the back end and its head, chunks, and end go back to the page", async () => {
    const chunk = new Uint8Array([104, 105]);
    const exchange = fakeExchange({ chunks: [chunk] });
    const { backend, calls } = fakeBackend(exchange);
    const endpoint = fakeEndpoint();

    serve(endpoint, backend);
    endpoint.receive({
      kind: "request",
      request: {
        seq: 3,
        method: "PUT",
        url: "/api/runs/r1/inputs/aln.fasta",
        headers: [
          { name: "content-type", value: "application/octet-stream" },
          { name: "x-a", value: "1" },
        ],
        body: new Uint8Array([1, 2]),
      },
    });
    await ended(endpoint, 3);

    expect([calls, endpoint.posted, exchange.freed()]).toStrictEqual([
      [
        [
          "PUT",
          "/api/runs/r1/inputs/aln.fasta",
          ["content-type", "x-a"],
          ["application/octet-stream", "1"],
          new Uint8Array([1, 2]),
        ],
      ],
      [
        {
          reply: { kind: "head", seq: 3, status: 200, headers: [{ name: "content-type", value: "text/plain" }] },
          transfer: [],
        },
        { reply: { kind: "chunk", seq: 3, data: chunk }, transfer: [chunk.buffer] },
        { reply: { kind: "end", seq: 3 }, transfer: [] },
      ],
      true,
    ]);
  });

  test("an abort from the page aborts the exchange and sends nothing more for it", async () => {
    const exchange = fakeExchange({ endless: true });
    const endpoint = fakeEndpoint();

    serve(endpoint, fakeBackend(exchange).backend);
    endpoint.receive(REQUEST);
    await pWaitFor(() => endpoint.posted.length === 1);
    endpoint.receive({ kind: "abort", seq: 7 });
    await pWaitFor(() => exchange.freed());

    expect([exchange.aborted(), endpoint.posted.map(({ reply }) => reply.kind)]).toStrictEqual([true, ["head"]]);
  });

  test("an error of the back end resets the exchange with its message", async () => {
    const exchange = fakeExchange({ chunkError: new Error("When reading the response body: gone") });
    const endpoint = fakeEndpoint();

    serve(endpoint, fakeBackend(exchange).backend);
    endpoint.receive(REQUEST);
    await ended(endpoint, 7);

    expect(endpoint.posted.map(({ reply }) => reply).at(-1)).toStrictEqual({
      kind: "reset",
      seq: 7,
      message: "When reading the response body: gone",
    });
  });

  test("a back end that throws when the request starts resets the exchange", async () => {
    const endpoint = fakeEndpoint();

    const backend: WasmBackend = {
      fetch: () => {
        throw new Error("a request has 2 header names and 1 header values");
      },
    };

    serve(endpoint, backend);
    endpoint.receive(REQUEST);
    await ended(endpoint, 7);

    expect(endpoint.posted).toStrictEqual([
      {
        reply: { kind: "reset", seq: 7, message: "a request has 2 header names and 1 header values" },
        transfer: [],
      },
    ]);
  });

  test("a trap of the module reports the trap once and resets every later request", async () => {
    const exchange = fakeExchange({ chunkError: new WebAssembly.RuntimeError("unreachable") });
    const endpoint = fakeEndpoint();
    const traps: WebAssembly.RuntimeError[] = [];

    serve(endpoint, fakeBackend(exchange).backend, (trap) => {
      traps.push(trap);
    });
    endpoint.receive(REQUEST);
    await ended(endpoint, 7);
    endpoint.receive({ ...REQUEST, request: { ...REQUEST.request, seq: 8 } });
    await ended(endpoint, 8);

    expect([traps.map(({ message }) => message), endpoint.posted.map(({ reply }) => reply).slice(1)]).toStrictEqual([
      ["unreachable"],
      [
        { kind: "reset", seq: 7, message: TRAPPED_MESSAGE },
        { kind: "reset", seq: 8, message: TRAPPED_MESSAGE },
      ],
    ]);
  });

  test("the job host runs after the head of a response is posted and before its body", async () => {
    const endpoint = fakeEndpoint();
    const steps: string[] = [];

    serve(endpoint, fakeBackend(fakeExchange({ chunks: [new Uint8Array([1])] })).backend, undefined, () => {
      steps.push(`job host after ${endpoint.posted.map(({ reply }) => reply.kind).join(",")}`);
    });
    endpoint.receive(REQUEST);
    await ended(endpoint, 7);

    expect(steps).toStrictEqual(["job host after head"]);
  });

  test("a message that is not a port message is refused", () => {
    const endpoint = fakeEndpoint();

    serve(endpoint, fakeBackend(fakeExchange({})).backend);

    expect(() => {
      endpoint.receive({ kind: "close", seq: 1 });
    }).toThrow("No matching discriminator");
  });
});

describe("serve response headers", () => {
  test("the header names and values of a response head become header pairs", () => {
    expect(responseHeaders({ status: 200, headerNames: ["a", "b"], headerValues: ["1", "2"] })).toStrictEqual([
      { name: "a", value: "1" },
      { name: "b", value: "2" },
    ]);
  });

  test("a response head with unequal names and values is refused", () => {
    expect(() => responseHeaders({ status: 200, headerNames: ["a", "b"], headerValues: ["1"] })).toThrow(
      "a response head has 2 header names and 1 header values",
    );
  });
});

interface FakeEndpoint extends ServeEndpoint {
  posted: Array<{ reply: PortReply; transfer: Transferable[] }>;
  receive(data: unknown): void;
}

interface FakeBackend {
  backend: WasmBackend;
  calls: Array<Parameters<WasmBackend["fetch"]>>;
}

interface FakeExchange extends WasmExchange {
  aborted(): boolean;
  freed(): boolean;
}

function serve(
  endpoint: ServeEndpoint,
  backend: WasmBackend,
  onTrap = (_trap: WebAssembly.RuntimeError) => undefined,
  afterHead = () => undefined,
): void {
  new ExchangeServer(endpoint, backend, onTrap, afterHead).listen();
}

function fakeEndpoint(): FakeEndpoint {
  let deliver: (data: unknown) => void = () => undefined;

  const endpoint: FakeEndpoint = {
    posted: [],
    post(reply, transfer) {
      endpoint.posted.push({ reply, transfer });
    },
    listen(listener) {
      deliver = listener;
    },
    receive(data) {
      deliver(data);
    },
  };

  return endpoint;
}

function fakeBackend(exchange: WasmExchange): FakeBackend {
  const calls: Array<Parameters<WasmBackend["fetch"]>> = [];

  return {
    calls,
    backend: {
      fetch: (...args) => {
        calls.push(args);

        return exchange;
      },
    },
  };
}

function fakeExchange({
  chunks = [],
  chunkError,
  endless = false,
}: {
  chunks?: Uint8Array[];
  chunkError?: Error;
  endless?: boolean;
}): FakeExchange {
  const head: WasmResponseHead = { status: 200, headerNames: ["content-type"], headerValues: ["text/plain"] };
  const remaining = [...chunks];
  const { promise: abortedRead, reject: rejectRead } = Promise.withResolvers<undefined>();
  let aborted = false;
  let freed = false;

  return {
    head: () => Promise.resolve(head),
    chunk: async () => {
      if (endless) {
        return abortedRead;
      }

      if (chunkError !== undefined) {
        throw chunkError;
      }

      return remaining.shift();
    },
    abort: () => {
      aborted = true;
      rejectRead(new Error("the exchange was aborted"));
    },
    free: () => {
      freed = true;
    },
    aborted: () => aborted,
    freed: () => freed,
  };
}

async function ended(endpoint: FakeEndpoint, seq: number): Promise<void> {
  await pWaitFor(() =>
    endpoint.posted.some(({ reply }) => reply.seq === seq && (reply.kind === "end" || reply.kind === "reset")),
  );
}
