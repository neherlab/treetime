import * as z from "zod";

const zSeq = z.int().nonnegative();

export const zBackendRequest = z.discriminatedUnion("kind", [
  z.object({ kind: z.literal("call"), seq: zSeq, request: z.string() }),
  z.object({ kind: z.literal("subscribe"), seq: zSeq, id: z.string(), from: z.int().nonnegative() }),
  z.object({ kind: z.literal("unsubscribe"), seq: zSeq }),
]);

export const zBackendReply = z.discriminatedUnion("kind", [
  z.object({ kind: z.literal("result"), seq: zSeq, json: z.string() }),
  z.object({ kind: z.literal("error"), seq: zSeq, error: z.string() }),
  z.object({ kind: z.literal("event"), seq: zSeq, json: z.string() }),
]);

export const zControlRequest = z.discriminatedUnion("kind", [
  z.object({ kind: z.literal("port") }),
  z.object({ kind: z.literal("save-file"), seq: zSeq, id: z.string(), path: z.string(), destination: z.string() }),
  z.object({ kind: z.literal("save-archive"), seq: zSeq, id: z.string(), destination: z.string() }),
]);

export const zControlReply = z.discriminatedUnion("kind", [
  z.object({ kind: z.literal("saved"), seq: zSeq }),
  z.object({ kind: z.literal("error"), seq: zSeq, error: z.string() }),
]);

export type BackendRequest = z.infer<typeof zBackendRequest>;

export type BackendReply = z.infer<typeof zBackendReply>;

export type ControlRequest = z.infer<typeof zControlRequest>;

export type SaveRequest = Exclude<ControlRequest, { kind: "port" }>;

export type ControlReply = z.infer<typeof zControlReply>;

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
