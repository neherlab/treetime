import { errorMessage, type ErrorCode } from "@neherlab/app-contracts";
import type { Backend, PortExchange, PortRequest, Subscription } from "@neherlab/app-napi";

import type { BackendRequest, ControlReply, FetchEndpoint, HostEndpoint, SaveRequest } from "./backend-protocol";

type Unsubscribe = Pick<Subscription, "unsubscribe">;

type Abort = Pick<PortExchange, "abort">;

export interface AddonBackend extends Pick<Backend, "call" | "saveRunFile" | "saveRunArchive"> {
  subscribe(...args: Parameters<Backend["subscribe"]>): Unsubscribe;
  fetch(...args: Parameters<Backend["fetch"]>): Abort;
}

const INVALID_REQUEST: ErrorCode = "invalid_request";

export function serveFetch(endpoint: FetchEndpoint, backend: Pick<AddonBackend, "fetch">): void {
  const exchanges = new Map<number, Abort>();

  const open = (request: PortRequest) => {
    const { seq } = request;

    try {
      let ended = false;

      // oxlint-disable-next-line treetime/require-io-timeout -- the addon runs the request in process; the renderer's signal ends it with an abort message
      const exchange = backend.fetch(request, (reply) => {
        if (reply.kind === "end" || reply.kind === "error") {
          ended = true;
          exchanges.delete(seq);
        }

        endpoint.post(reply);
      });

      if (!ended) {
        exchanges.set(seq, exchange);
      }
    } catch (error: unknown) {
      endpoint.post({ kind: "error", seq, error: { code: INVALID_REQUEST, message: errorMessage(error), causes: [] } });
    }
  };

  endpoint.listen((message) => {
    if (message.kind === "request") {
      open(message.request);
    } else {
      exchanges.get(message.seq)?.abort();
      exchanges.delete(message.seq);
    }
  });

  endpoint.onClose(() => {
    for (const exchange of exchanges.values()) {
      exchange.abort();
    }

    exchanges.clear();
  });
}

export function serveBackend(endpoint: HostEndpoint, backend: AddonBackend): void {
  const subscriptions = new Map<number, Unsubscribe>();

  const answer = async (request: BackendRequest): Promise<void> => {
    switch (request.kind) {
      case "call":
        endpoint.post({ kind: "result", seq: request.seq, json: await backend.call(request.request) });
        break;
      case "subscribe":
        subscriptions.set(
          request.seq,
          backend.subscribe(request.id, request.from, (err, json) => {
            if (err === null) {
              endpoint.post({ kind: "event", seq: request.seq, json });
            }
          }),
        );
        break;
      case "unsubscribe":
        subscriptions.get(request.seq)?.unsubscribe();
        subscriptions.delete(request.seq);
        break;
    }
  };

  const answerOrReport = async (request: BackendRequest) => {
    try {
      await answer(request);
    } catch (error: unknown) {
      endpoint.post({ kind: "error", seq: request.seq, error: errorMessage(error) });
    }
  };

  endpoint.listen((request) => void answerOrReport(request));

  endpoint.onClose(() => {
    for (const subscription of subscriptions.values()) {
      subscription.unsubscribe();
    }

    subscriptions.clear();
  });
}

export async function saveRunFiles(backend: AddonBackend, request: SaveRequest): Promise<ControlReply> {
  try {
    await (request.kind === "save-file"
      ? backend.saveRunFile(request.id, request.path, request.destination)
      : backend.saveRunArchive(request.id, request.destination));

    return { kind: "saved", seq: request.seq };
  } catch (error: unknown) {
    return { kind: "error", seq: request.seq, error: errorMessage(error) };
  }
}
