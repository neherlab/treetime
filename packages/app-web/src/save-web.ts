import { type ApiClient, runsArchive, runsFile } from "@neherlab/app-contracts/client";
import type { SaveActions } from "@neherlab/app-ui";

export function createWebSaveActions(
  client: ApiClient,
  saveBlob: (blob: Blob, name: string) => void = downloadBlob,
): SaveActions {
  return {
    async saveRunFile(id, path, name) {
      const { data } = await runsFile({ client, path: { id }, query: { path }, parseAs: "blob", throwOnError: true });
      saveBlob(data, name);

      return true;
    },
    async saveRunArchive(id, name) {
      const { data } = await runsArchive({ client, path: { id }, parseAs: "blob", throwOnError: true });
      saveBlob(data, name);

      return true;
    },
  };
}

export function downloadBlob(blob: Blob, name: string): void {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");

  anchor.href = url;
  anchor.download = name;
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  URL.revokeObjectURL(url);
}
