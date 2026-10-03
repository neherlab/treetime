import {
  zErrorResponse,
  zPickedFiles,
  zPickedFolder,
  type ErrorResponse,
  type LocalFiles,
  type PickFilesRequest,
  type PickFolderRequest,
  type WorkspaceShell,
} from "@neherlab/app-contracts";

import { BACKEND_PORT_CHANNEL } from "./channels";
import type { FetchConnection, FetchPort } from "./port-fetch";
import type { SaveReply, SaveRunArchiveDialog, SaveRunFileDialog } from "./shell-protocol";

interface WindowMessage {
  source: unknown;
  data: unknown;
  ports: readonly [FetchPort?];
}

export interface WindowLike {
  addEventListener(type: "message", listener: (event: WindowMessage) => void): void;
}

export interface DesktopShell {
  connectBackend(): void;
  onBackendStopped(listener: (reason: string, restarts: boolean) => void): void;
  pickFiles(request: PickFilesRequest): Promise<unknown>;
  pickFolder(request: PickFolderRequest): Promise<unknown>;
  restartBackend(): Promise<unknown>;
  pathForFile(file: File): string;
  saveRunFile(request: SaveRunFileDialog): Promise<SaveReply>;
  saveRunArchive(request: SaveRunArchiveDialog): Promise<SaveReply>;
}

export class SaveError extends Error {
  readonly response: ErrorResponse;

  constructor(response: ErrorResponse) {
    super([response.message, ...response.causes].join(": "));
    this.name = "SaveError";
    this.response = response;
  }
}

export function createDesktopSaveActions(shell: DesktopShell) {
  return {
    saveRunFile: async (id: string, path: string, name: string) => saved(await shell.saveRunFile({ id, path, name })),
    saveRunArchive: async (id: string, name: string) => saved(await shell.saveRunArchive({ id, name })),
  };
}

export function createLocalFiles(shell: DesktopShell): LocalFiles {
  return {
    async pickFiles(request) {
      return zPickedFiles.parse(await shell.pickFiles(request));
    },
    pathForFile: (file) => shell.pathForFile(file),
  };
}

export function createWorkspaceShell(shell: DesktopShell): WorkspaceShell {
  return {
    async pickFolder(request) {
      return zPickedFolder.parse(await shell.pickFolder(request));
    },
    async restartBackend() {
      await shell.restartBackend();
    },
  };
}

export function windowFetchConnection(target: WindowLike, shell: DesktopShell): FetchConnection {
  const portListeners: Array<(port: FetchPort) => void> = [];
  let current: FetchPort | undefined;
  let requested = false;

  shell.onBackendStopped(() => {
    current = undefined;
  });

  target.addEventListener("message", (event) => {
    const [port] = event.ports;

    if (event.source !== target || !isPortMessage(event.data) || port === undefined) {
      return;
    }

    current = port;
    portListeners.forEach((listener) => {
      listener(port);
    });
  });

  return {
    onPort(listener) {
      portListeners.push(listener);

      if (current !== undefined) {
        listener(current);
      } else if (!requested) {
        requested = true;
        shell.connectBackend();
      }
    },
    onStopped(listener) {
      shell.onBackendStopped(listener);
    },
  };
}

function isPortMessage(data: unknown): boolean {
  return typeof data === "object" && data !== null && "channel" in data && data.channel === BACKEND_PORT_CHANNEL;
}

function saved(reply: SaveReply): boolean {
  if ("error" in reply) {
    throw new SaveError(errorResponse(reply.error));
  }

  return reply.saved;
}

function errorResponse(text: string): ErrorResponse {
  const parsed = zErrorResponse.safeParse(parseJson(text));

  return parsed.success ? parsed.data : { code: "internal_error", message: text, causes: [] };
}

function parseJson(text: string): unknown {
  try {
    const value: unknown = JSON.parse(text);

    return value;
  } catch {
    return undefined;
  }
}
