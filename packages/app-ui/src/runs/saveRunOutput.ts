import { runsArchive, runsFile, type ApiClient } from "@neherlab/app-contracts/client";

import { requestUrl } from "../api/keys";
import type { Host, SaveRunReply } from "../host";

export interface RunOutput {
  id: string;
  path?: string;
  name: string;
}

export interface SaveTarget {
  host: Pick<Host, "saveRun"> | null;
  client: ApiClient;
  download: (url: string, name: string) => void;
}

export async function saveRunOutput({ host, client, download }: SaveTarget, output: RunOutput): Promise<SaveRunReply> {
  if (host !== null) {
    return host.saveRun(output);
  }

  const { id, path, name } = output;

  const url =
    path === undefined
      ? requestUrl(client, (context) => runsArchive({ ...context, path: { id } }))
      : requestUrl(client, (context) => runsFile({ ...context, path: { id }, query: { path } }));

  download(url, name);

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
