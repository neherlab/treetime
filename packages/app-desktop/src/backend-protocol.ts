import * as z from "zod";

const zSeq = z.int().nonnegative();

export const zBackendRequest = z.discriminatedUnion("kind", [
  z.object({ kind: z.literal("call"), seq: zSeq, request: z.string() }),
  z.object({ kind: z.literal("subscribe"), seq: zSeq, id: z.string(), from: z.int().nonnegative() }),
  z.object({ kind: z.literal("unsubscribe"), seq: zSeq }),
  z.object({ kind: z.literal("read-file"), seq: zSeq, id: z.string(), path: z.string() }),
  z.object({ kind: z.literal("archive"), seq: zSeq, id: z.string() }),
]);

export const zBackendReply = z.discriminatedUnion("kind", [
  z.object({ kind: z.literal("result"), seq: zSeq, json: z.string() }),
  z.object({ kind: z.literal("error"), seq: zSeq, error: z.string() }),
  z.object({ kind: z.literal("event"), seq: zSeq, json: z.string() }),
  z.object({ kind: z.literal("chunk"), seq: zSeq, bytes: z.instanceof(ArrayBuffer) }),
  z.object({ kind: z.literal("end"), seq: zSeq }),
]);

export type BackendRequest = z.infer<typeof zBackendRequest>;

export type BackendReply = z.infer<typeof zBackendReply>;

export interface MessageEndpoint<Incoming, Outgoing> {
  post(message: Outgoing): void;
  listen(listener: (message: Incoming) => void): void;
  onClose(listener: () => void): void;
}

export type HostEndpoint = MessageEndpoint<BackendRequest, BackendReply>;

export type ClientEndpoint = MessageEndpoint<BackendReply, BackendRequest>;

export interface PortLike {
  postMessage(message: BackendRequest | BackendReply): void;
  addEventListener(type: "message" | "close", listener: (event: { data: unknown }) => void): void;
  start(): void;
}

export function portEndpoint<Incoming, Outgoing extends BackendRequest | BackendReply>(
  port: PortLike,
  schema: z.ZodType<Incoming>,
): MessageEndpoint<Incoming, Outgoing> {
  return {
    post: (message) => {
      port.postMessage(message);
    },
    listen: (listener) => {
      port.addEventListener("message", (event) => {
        const message = schema.safeParse(event.data);

        if (message.success) {
          listener(message.data);
        } else {
          console.warn("[TreeTime] ignored a malformed message", message.error.message);
        }
      });
      port.start();
    },
    onClose: (listener) => {
      port.addEventListener("close", () => {
        listener();
      });
    },
  };
}
