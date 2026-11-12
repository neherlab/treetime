import { runsFile, runsPackage, type ApiClient } from "@neherlab/app-contracts/client";

import type { Host, SaveRunReply } from "../host";

const OBJECT_URL_LIFETIME_MS = 40_000;

const PICKER_CANCELED = "AbortError";

export interface RunOutput {
  id: string;
  path?: string;
  name: string;
}

export interface SaveTarget {
  host: Pick<Host, "saveRun"> | null;
  client: ApiClient;
  pickFile: SaveFilePicker | null;
  download: (url: string, name: string) => void;
  objectUrls: ObjectUrls;
  setTimer: (callback: () => void, delayMs: number) => void;
}

export type SaveFilePicker = (name: string) => Promise<WritableStream<Uint8Array> | undefined>;

export interface ObjectUrls {
  createObjectURL(blob: Blob): string;
  revokeObjectURL(url: string): void;
}

export async function saveRunOutput(
  { host, client, pickFile, download, objectUrls, setTimer }: SaveTarget,
  output: RunOutput,
): Promise<SaveRunReply> {
  if (host !== null) {
    return host.saveRun(output);
  }

  if (pickFile !== null) {
    return streamToFile(client, pickFile, output);
  }

  const { data } = await fetchOutput(client, output, "blob");
  const url = objectUrls.createObjectURL(data);

  try {
    download(url, output.name);
  } finally {
    setTimer(() => {
      objectUrls.revokeObjectURL(url);
    }, OBJECT_URL_LIFETIME_MS);
  }

  return { kind: "saved" };
}

export function browserFilePicker(scope: SavePickerScope = window): SaveFilePicker | null {
  const showSaveFilePicker = scope.showSaveFilePicker;

  if (!scope.isSecureContext || showSaveFilePicker === undefined) {
    return null;
  }

  return async (name) => {
    try {
      const handle = await showSaveFilePicker({ suggestedName: name });

      return await handle.createWritable();
    } catch (error: unknown) {
      if (error instanceof DOMException && error.name === PICKER_CANCELED) {
        return undefined;
      }

      throw error;
    }
  };
}

export function downloadUrl(url: string, name: string): void {
  const anchor = document.createElement("a");

  anchor.href = url;
  anchor.download = name;
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
}

export function setBrowserTimer(callback: () => void, delayMs: number): void {
  globalThis.setTimeout(callback, delayMs);
}

export function saveText(text: string, name: string, type: string): void {
  const url = URL.createObjectURL(new Blob([text], { type }));

  try {
    downloadUrl(url, name);
  } finally {
    setBrowserTimer(() => URL.revokeObjectURL(url), OBJECT_URL_LIFETIME_MS);
  }
}

async function streamToFile(client: ApiClient, pickFile: SaveFilePicker, output: RunOutput): Promise<SaveRunReply> {
  const sink = await pickFile(output.name);

  if (sink === undefined) {
    return { kind: "canceled" };
  }

  let body: ReadableStream<Uint8Array> | null;

  try {
    ({
      response: { body },
    } = await fetchOutput(client, output, "stream"));
  } catch (error: unknown) {
    await sink.abort(error);
    throw error;
  }

  if (body === null) {
    await sink.close();
  } else {
    await body.pipeTo(sink);
  }

  return { kind: "saved" };
}

async function fetchOutput(client: ApiClient, { id, path }: RunOutput, parseAs: "blob" | "stream") {
  return path === undefined
    ? runsPackage({ client, path: { id }, parseAs, throwOnError: true })
    : runsFile({ client, path: { id }, query: { path }, parseAs, throwOnError: true });
}

interface SavePickerScope {
  isSecureContext: boolean;
  showSaveFilePicker?: (options: { suggestedName: string }) => Promise<FileSystemFileHandle>;
}
