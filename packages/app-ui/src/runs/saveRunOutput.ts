import { runsFile, runsPackage, type ApiClient } from "@neherlab/app-contracts/client";

import type { Host, SaveRunReply } from "../host";

const OBJECT_URL_LIFETIME_MS = 40_000;

export interface RunOutput {
  id: string;
  path?: string;
  name: string;
}

export interface SaveTarget {
  host: Pick<Host, "saveRun"> | null;
  client: ApiClient;
  download: (url: string, name: string) => void;
  objectUrls: ObjectUrls;
  setTimer: (callback: () => void, delayMs: number) => void;
}

export interface ObjectUrls {
  createObjectURL(blob: Blob): string;
  revokeObjectURL(url: string): void;
}

export async function saveRunOutput(
  { host, client, download, objectUrls, setTimer }: SaveTarget,
  output: RunOutput,
): Promise<SaveRunReply> {
  if (host !== null) {
    return host.saveRun(output);
  }

  const { id, path, name } = output;

  const { data } =
    path === undefined
      ? await runsPackage({ client, path: { id }, parseAs: "blob", throwOnError: true })
      : await runsFile({ client, path: { id }, query: { path }, parseAs: "blob", throwOnError: true });

  const url = objectUrls.createObjectURL(data);

  try {
    download(url, name);
  } finally {
    setTimer(() => {
      objectUrls.revokeObjectURL(url);
    }, OBJECT_URL_LIFETIME_MS);
  }

  return { kind: "saved" };
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
