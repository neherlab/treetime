import { errorMessage, type AppCommand } from "@neherlab/app-contracts";
import { runsCreate, runsReleaseInput, runsUploadInput } from "@neherlab/app-contracts/client";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useCallback, useMemo } from "react";

import { useApiContext } from "../api/context";
import { useDraftStore } from "../store/draft";
import { type StoredFile, storeInOrder } from "./fileSelection";
import { pendingUploadRun } from "./pendingUpload";

export interface RunUploads {
  uploadRun: (command: AppCommand) => Promise<string>;
  uploadFiles: (runId: string, files: readonly File[], replaces: readonly string[]) => Promise<StoredFile[]>;
  release: (runId: string, paths: readonly string[]) => void;
  releaseFromPendingRun: (runId: string, paths: readonly string[]) => void;
}

export function useRunUploads(): RunUploads {
  const { client } = useApiContext();
  const toasts = useToastManager();

  const releaseOne = useCallback(
    async (runId: string, path: string) => {
      try {
        await runsReleaseInput({ client, path: { id: runId }, query: { path }, throwOnError: true });
      } catch (error) {
        toasts.add({ title: "An earlier upload cannot be removed", description: errorMessage(error) });
      }
    },
    [client, toasts],
  );

  const release = useCallback(
    (runId: string, paths: readonly string[]) => {
      for (const path of paths) {
        void releaseOne(runId, path);
      }
    },
    [releaseOne],
  );

  const releaseFromPendingRun = useCallback(
    (runId: string, paths: readonly string[]) => {
      if (useDraftStore.getState().draft.upload_run_id === runId) {
        release(runId, paths);
      }
    },
    [release],
  );

  const uploadRun = useCallback(
    async (command: AppCommand): Promise<string> => {
      const pending = useDraftStore.getState().uploadRunCreation;

      if (pending !== undefined) {
        return pending;
      }

      const existing = await pendingUploadRun(client, useDraftStore.getState().draft.upload_run_id);

      if (existing !== undefined) {
        return existing.id;
      }

      const started = useDraftStore.getState().uploadRunCreation;

      if (started !== undefined) {
        return started;
      }

      const creation: Promise<string> = runsCreate({
        client,
        body: { command, config: {}, defer_start: true },
        throwOnError: true,
      })
        .then(({ data: record }) => {
          useDraftStore.getState().settleUploadRun(creation, record.id);

          return record.id;
        })
        .finally(() => useDraftStore.getState().settleUploadRun(creation, undefined));

      useDraftStore.getState().beginUploadRun(creation);

      return creation;
    },
    [client],
  );

  const uploadFiles = useCallback(
    async (runId: string, files: readonly File[], replaces: readonly string[]) =>
      storeInOrder(
        files,
        async (file) => {
          const { data: uploaded } = await runsUploadInput({
            client,
            path: { id: runId, name: file.name },
            query: { replaces: [...replaces] },
            body: file,
            throwOnError: true,
          });

          return { name: file.name, path: uploaded.path, size: uploaded.size };
        },
        (partial) =>
          releaseFromPendingRun(
            runId,
            partial.map((file) => file.path),
          ),
      ),
    [client, releaseFromPendingRun],
  );

  return useMemo(
    () => ({ uploadRun, uploadFiles, release, releaseFromPendingRun }),
    [release, releaseFromPendingRun, uploadFiles, uploadRun],
  );
}
