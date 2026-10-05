import {
  BACKEND_PORT_CHANNEL,
  BACKEND_PORT_REQUEST_CHANNEL,
  hostChannelName,
  hostChannelSchemas,
  hostEventSchema,
  type HostChannel,
  type HostEvent,
  type HostEventPayload,
  type HostHandlerReply,
  type HostHandlerRequest,
} from "@neherlab/app-ui/host";

import { isTrustedSender } from "./security";

export interface SenderEvent<Sender> {
  sender: Sender;
  senderFrame: { url: string; parent: unknown } | null;
}

export interface IpcMainLike<Sender> {
  // oxlint-disable-next-line anti-slop/no-unknown-returns -- Electron IPC requests and replies are untyped; handle() parses both with the channel schema
  handle(channel: string, listener: (event: SenderEvent<Sender>, ...args: unknown[]) => unknown): void;
  on(channel: string, listener: (event: SenderEvent<Sender>, ...args: unknown[]) => void): void;
}

export interface WebContentsLike<Port> {
  send(channel: string, ...args: unknown[]): void;
  postMessage(channel: string, message: null, transfer: Port[]): void;
}

export function handle<Sender, C extends HostChannel>(
  ipc: IpcMainLike<Sender>,
  appUrl: string,
  channel: C,
  handler: (sender: Sender, request: HostHandlerRequest<C>) => Promise<HostHandlerReply<C>> | HostHandlerReply<C>,
): void {
  const name = hostChannelName(channel);
  const schemas = hostChannelSchemas(channel);

  ipc.handle(name, async (event, ...args) => {
    if (!isTrustedSender(event.senderFrame, appUrl)) {
      throw new Error(`${name} refused a message from a frame outside the application`);
    }

    return schemas.reply.parse(await handler(event.sender, schemas.request.parse(args[0])));
  });
}

export function listenForPortRequests<Sender>(
  ipc: IpcMainLike<Sender>,
  appUrl: string,
  listener: (sender: Sender) => void,
): void {
  ipc.on(BACKEND_PORT_REQUEST_CHANNEL, (event) => {
    if (isTrustedSender(event.senderFrame, appUrl)) {
      listener(event.sender);
    }
  });
}

export function emit<E extends HostEvent, Port>(
  contents: WebContentsLike<Port>,
  event: E,
  payload: HostEventPayload<E>,
): void {
  contents.send(hostChannelName(event), hostEventSchema(event).parse(payload));
}

export function sendPort<Port>(contents: WebContentsLike<Port>, port: Port): void {
  contents.postMessage(BACKEND_PORT_CHANNEL, null, [port]);
}
