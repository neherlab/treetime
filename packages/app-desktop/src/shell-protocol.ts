import * as z from "zod";

import { BACKEND_PORT_CHANNEL, BACKEND_STOPPED_CHANNEL } from "./channels";

export const zShellMessage = z.discriminatedUnion("channel", [
  z.object({ channel: z.literal(BACKEND_PORT_CHANNEL) }),
  z.object({ channel: z.literal(BACKEND_STOPPED_CHANNEL), reason: z.string(), restarts: z.boolean() }),
]);

export type ShellMessage = z.infer<typeof zShellMessage>;

export const zSaveRunFileRequest = z.object({ id: z.string(), path: z.string(), name: z.string() });

export type SaveRunFileRequest = z.infer<typeof zSaveRunFileRequest>;

export const zSaveRunArchiveRequest = z.object({ id: z.string(), name: z.string() });

export type SaveRunArchiveRequest = z.infer<typeof zSaveRunArchiveRequest>;

export const zSaveReply = z.union([z.object({ saved: z.boolean() }), z.object({ error: z.string() })]);

export type SaveReply = z.infer<typeof zSaveReply>;
