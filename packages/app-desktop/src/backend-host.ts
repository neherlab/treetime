import type { Backend, Subscription } from "@neherlab/app-napi";

import type { BackendRequest, HostEndpoint } from "./backend-protocol";

type Unsubscribe = Pick<Subscription, "unsubscribe">;

export interface AddonBackend extends Pick<Backend, "call" | "readFile" | "archive"> {
  subscribe(...args: Parameters<Backend["subscribe"]>): Unsubscribe;
}

const CHUNK_SIZE = 1 << 20;

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
      case "read-file":
        sendBytes(endpoint, request.seq, backend.readFile(request.id, request.path));
        break;
      case "archive":
        sendBytes(endpoint, request.seq, backend.archive(request.id));
        break;
    }
  };

  endpoint.listen((request) => {
    answer(request).catch((error: Error) => {
      endpoint.post({ kind: "error", seq: request.seq, error: error.message });
    });
  });

  endpoint.onClose(() => {
    for (const subscription of subscriptions.values()) {
      subscription.unsubscribe();
    }

    subscriptions.clear();
  });
}

function sendBytes(endpoint: HostEndpoint, seq: number, bytes: Uint8Array): void {
  for (let offset = 0; offset < bytes.byteLength; offset += CHUNK_SIZE) {
    const chunk = new Uint8Array(bytes.subarray(offset, offset + CHUNK_SIZE));
    endpoint.post({ kind: "chunk", seq, bytes: chunk.buffer });
  }

  endpoint.post({ kind: "end", seq });
}
