import { zErrorResponse, zRunsSavePath, zSaveRunRequest, zUiTheme } from "@neherlab/app-contracts";
import * as z from "zod";

export const HOST_CHANNEL_PREFIX = "treetime:";

export const HOST_CHANNELS = {
  "pick-files": {
    request: z.strictObject({ title: z.string(), extensions: z.array(z.string()), multiple: z.boolean() }),
    reply: z.array(z.string()),
  },
  "pick-folder": {
    request: z.strictObject({ title: z.string() }),
    reply: z.string().optional(),
  },
  "save-run": {
    request: z.strictObject({ id: zRunsSavePath.shape.id, path: zSaveRunRequest.shape.path, name: z.string() }),
    reply: z.discriminatedUnion("kind", [
      z.strictObject({ kind: z.literal("saved") }),
      z.strictObject({ kind: z.literal("canceled") }),
      z.strictObject({ kind: z.literal("error"), error: zErrorResponse }),
    ]),
  },
  "restart-backend": {
    request: z.undefined(),
    reply: z.undefined(),
  },
  "set-native-theme": {
    request: zUiTheme,
    reply: z.undefined(),
  },
} as const;

export const HOST_EVENTS = {
  "backend-stopped": z.strictObject({
    reason: z.string(),
    restarts: z.boolean(),
    error: zErrorResponse.optional(),
  }),
} as const;

export const BACKEND_PORT_REQUEST_CHANNEL = `${HOST_CHANNEL_PREFIX}backend-port-request`;

export const BACKEND_PORT_CHANNEL = `${HOST_CHANNEL_PREFIX}backend-port`;

export type HostChannel = keyof typeof HOST_CHANNELS;

export type HostRequest<C extends HostChannel> = z.input<(typeof HOST_CHANNELS)[C]["request"]>;

export type HostReply<C extends HostChannel> = z.output<(typeof HOST_CHANNELS)[C]["reply"]>;

export type HostEvent = keyof typeof HOST_EVENTS;

export type HostEventPayload<E extends HostEvent> = z.output<(typeof HOST_EVENTS)[E]>;

export type BackendStopped = HostEventPayload<"backend-stopped">;

export type SaveRunReply = HostReply<"save-run">;

export type HostHandlerRequest<C extends HostChannel> = z.output<(typeof HOST_CHANNELS)[C]["request"]>;

export type HostHandlerReply<C extends HostChannel> = z.input<(typeof HOST_CHANNELS)[C]["reply"]>;

export interface HostChannelSchemas<C extends HostChannel> {
  request: z.ZodType<HostHandlerRequest<C>, HostRequest<C>>;
  reply: z.ZodType<HostReply<C>, HostHandlerReply<C>>;
}

export interface Host {
  pickFiles(request: HostRequest<"pick-files">): Promise<HostReply<"pick-files">>;
  pickFolder(request: HostRequest<"pick-folder">): Promise<HostReply<"pick-folder">>;
  saveRun(request: HostRequest<"save-run">): Promise<SaveRunReply>;
  restartBackend(): Promise<HostReply<"restart-backend">>;
  setNativeTheme(theme: HostRequest<"set-native-theme">): Promise<HostReply<"set-native-theme">>;
  pathForFile(file: File): string;
  connectBackend(): void;
  onBackendStopped(listener: (stop: BackendStopped) => void): void;
}

export function hostChannelName(channel: HostChannel | HostEvent): string {
  return `${HOST_CHANNEL_PREFIX}${channel}`;
}

export function hostChannelSchemas<C extends HostChannel>(channel: C): HostChannelSchemas<C> {
  return HOST_CHANNEL_TABLE[channel];
}

export function hostEventSchema<E extends HostEvent>(event: E): z.ZodType<HostEventPayload<E>> {
  return HOST_EVENT_TABLE[event];
}

type HostChannelTable = { [K in HostChannel]: HostChannelSchemas<K> };

type HostEventTable = { [K in HostEvent]: z.ZodType<HostEventPayload<K>> };

const HOST_CHANNEL_TABLE: HostChannelTable = HOST_CHANNELS;

const HOST_EVENT_TABLE: HostEventTable = HOST_EVENTS;
