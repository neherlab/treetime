import { create } from "zustand";

import type { BackendStopped } from "../host";

export type StopAction = "dismiss" | "reload-page" | "restart-back-end";

export interface BackendStopNotice {
  message: string;
  restarts: boolean;
}

export function stopAction(restarts: boolean, desktop: boolean): StopAction {
  if (restarts) {
    return "dismiss";
  }

  return desktop ? "restart-back-end" : "reload-page";
}

export function stopEffect(desktop: boolean): string {
  return desktop
    ? "Runs that were computing stopped. Finished runs stay in the runs folder."
    : "Every run of this tab is lost. Download run packages to keep runs.";
}

interface AppNoticesState {
  stop: BackendStopNotice | undefined;
  stops: number;
  chunkFailed: boolean;
}

interface AppNoticesActions {
  stopped: (stop: BackendStopped) => void;
  dismissStop: () => void;
  failChunk: () => void;
}

export const useAppNoticesStore = create<AppNoticesState & AppNoticesActions>()((set) => ({
  stop: undefined,
  stops: 0,
  chunkFailed: false,
  stopped: (stop) => {
    if (stop.requested) {
      return;
    }

    set(({ stops }) => ({ stop: { message: stopMessage(stop), restarts: stop.restarts }, stops: stops + 1 }));
  },
  dismissStop: () => set({ stop: undefined }),
  failChunk: () => set({ chunkFailed: true }),
}));

export function publishBackendStop(stop: BackendStopped): void {
  useAppNoticesStore.getState().stopped(stop);
}

function stopMessage({ reason, error }: BackendStopped): string {
  const text = error === undefined ? reason : [error.message, ...error.causes].join(": ");

  return text.split("\n", 1)[0] ?? text;
}
