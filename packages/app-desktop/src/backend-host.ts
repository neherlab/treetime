import type { Backend, Subscription } from "@neherlab/app-napi";
import * as z from "zod";

import type { BackendRequest, ControlReply, HostEndpoint, SaveRequest } from "./backend-protocol";

type Unsubscribe = Pick<Subscription, "unsubscribe">;

export interface AddonBackend extends Pick<Backend, "call" | "resolveRunFile" | "saveRunFile" | "saveRunArchive"> {
  subscribe(...args: Parameters<Backend["subscribe"]>): Unsubscribe;
}

export type ReadChunks = (path: string) => AsyncIterable<Uint8Array> | Iterable<Uint8Array>;

export function serveBackend(endpoint: HostEndpoint, backend: AddonBackend, readChunks: ReadChunks): void {
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
      case "read-file": {
        const path = z.string().parse(JSON.parse(await backend.resolveRunFile(request.id, request.path)));

        for await (const chunk of readChunks(path)) {
          endpoint.post({ kind: "chunk", seq: request.seq, bytes: new Uint8Array(chunk).buffer });
        }

        endpoint.post({ kind: "end", seq: request.seq });
        break;
      }
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

export async function saveRunFiles(backend: AddonBackend, request: SaveRequest): Promise<ControlReply> {
  try {
    await (request.kind === "save-file"
      ? backend.saveRunFile(request.id, request.path, request.destination)
      : backend.saveRunArchive(request.id, request.destination));

    return { kind: "saved", seq: request.seq };
  } catch (error: unknown) {
    return { kind: "error", seq: request.seq, error: error instanceof Error ? error.message : String(error) };
  }
}
