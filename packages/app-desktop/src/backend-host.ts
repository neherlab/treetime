import { errorMessage, type ErrorCode } from "@neherlab/app-contracts";
import type { Backend, PortExchange, PortRequest } from "@neherlab/app-napi";

import type { ControlReply, FetchEndpoint, SaveRequest } from "./backend-protocol";

type Abort = Pick<PortExchange, "abort">;

export interface AddonBackend extends Pick<Backend, "saveRunFile" | "saveRunArchive"> {
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

export async function saveRunFiles(backend: AddonBackend, request: SaveRequest): Promise<ControlReply> {
  try {
    await (request.kind === "save-file"
      ? backend.saveRunFile(request.request)
      : backend.saveRunArchive(request.request));

    return { kind: "saved", seq: request.seq };
  } catch (error: unknown) {
    return { kind: "error", seq: request.seq, error: errorMessage(error) };
  }
}
