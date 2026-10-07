import { errorMessage } from "@neherlab/app-contracts";
import type { Backend, PortExchange, PortRequest } from "@neherlab/app-napi";

import type { FetchEndpoint, PortScope } from "./backend-protocol";

type Abort = Pick<PortExchange, "abort">;

export interface AddonBackend {
  fetch(...args: Parameters<Backend["fetch"]>): Abort;
  rejectRequest: Backend["rejectRequest"];
}

export function serveFetch(endpoint: FetchEndpoint, backend: AddonBackend, scope: PortScope): void {
  const exchanges = new Map<number, Abort>();

  const open = (request: PortRequest) => {
    const { seq } = request;

    try {
      let ended = false;

      // oxlint-disable-next-line custom/require-io-timeout -- the addon runs the request in process; the renderer's signal ends it with an abort message
      const exchange = backend.fetch(request, scope, (reply) => {
        if (reply.kind === "end" || reply.kind === "reset") {
          ended = true;
          exchanges.delete(seq);
        }

        endpoint.post(reply);
      });

      if (!ended) {
        exchanges.set(seq, exchange);
      }
    } catch (error: unknown) {
      backend.rejectRequest(seq, errorMessage(error)).forEach((reply) => {
        endpoint.post(reply);
      });
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
