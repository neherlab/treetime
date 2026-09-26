import type { SaveRunArchiveRequest, SaveRunFileRequest } from "@neherlab/app-napi";
import * as z from "zod";

import { BACKEND_PORT_CHANNEL, BACKEND_STOPPED_CHANNEL } from "./channels";

export const zShellMessage = z.discriminatedUnion("channel", [
  z.object({ channel: z.literal(BACKEND_PORT_CHANNEL) }),
  z.object({ channel: z.literal(BACKEND_STOPPED_CHANNEL), reason: z.string(), restarts: z.boolean() }),
]);

export type ShellMessage = z.infer<typeof zShellMessage>;

export type SaveRunFileDialog = Omit<SaveRunFileRequest, "destination"> & { name: string };

export type SaveRunArchiveDialog = Omit<SaveRunArchiveRequest, "destination"> & { name: string };

export type SaveReply = { saved: boolean } | { error: string };
