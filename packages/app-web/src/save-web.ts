import { type ApiClient, runsArchive, runsFile } from "@neherlab/app-contracts/client";
import { requestUrl, type SaveActions } from "@neherlab/app-ui";

export function createWebSaveActions(
  client: ApiClient,
  saveUrl: (url: string, name: string) => void = downloadUrl,
): SaveActions {
  return {
    saveRunFile(id, path, name) {
      saveUrl(
        requestUrl(client, (context) => runsFile({ ...context, path: { id }, query: { path } })),
        name,
      );

      return Promise.resolve(true);
    },
    saveRunArchive(id, name) {
      saveUrl(
        requestUrl(client, (context) => runsArchive({ ...context, path: { id } })),
        name,
      );

      return Promise.resolve(true);
    },
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
