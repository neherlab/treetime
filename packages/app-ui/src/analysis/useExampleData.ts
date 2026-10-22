import type { AppCommand, IndexedConfig, IndexedDataset } from "@neherlab/app-contracts";
import { examples } from "@neherlab/app-contracts/client";
import { useQuery } from "@tanstack/react-query";
import { useNavigate } from "@tanstack/react-router";
import { useCallback, useMemo } from "react";

import { useApi } from "../api/hooks";
import { useHost } from "../host-context";
import { inputSettingKeys } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { applyDraftInputs } from "./draftInputs";
import {
  type ConfigInput,
  configInputs,
  DATA_INDEX_FILE,
  indexedDatasetInputs,
  type ExampleIndex,
  exampleDataUrl,
  fetchExampleFile,
  fetchExampleIndex,
  localExampleIndex,
  localExamplePath,
  localInputAssignment,
  uploadedConfigInputs,
  webExampleIndex,
} from "./exampleData";
import { draftUploadPaths, type StoredFile, storedSource, storedValue } from "./fileSelection";
import { claimInputs } from "./inputClaims";
import { type ConfigLoadResult, SUPERSEDED, useConfigLoader } from "./useConfigLoader";
import { useRunUploads } from "./useRunUploads";

export function useExampleIndex(): ExampleIndexQuery {
  const host = useHost();
  const base = exampleDataBase();

  const web = useQuery({
    queryKey: ["example-data", exampleDataUrl(base, DATA_INDEX_FILE)],
    queryFn: async ({ signal }) => fetchExampleIndex(base, signal),
    staleTime: Infinity,
    enabled: host === null,
  });

  const local = useApi((context) => examples(context), { staleTime: Infinity, enabled: host !== null });

  return useMemo(() => {
    if (host === null) {
      return { data: web.data === undefined ? undefined : webExampleIndex(web.data), error: web.error };
    }

    return { data: local.data === undefined ? undefined : localExampleIndex(local.data), error: local.error };
  }, [host, local.data, local.error, web.data, web.error]);
}

export function useExampleConfigLoader() {
  const loadConfig = useConfigLoader();
  const { uploadRun, uploadFiles, releaseFromPendingRun } = useRunUploads();

  return useCallback(
    async (examples: ExampleIndex, config: IndexedConfig): Promise<ConfigLoadResult> => {
      if (examples.origin === "local") {
        return loadConfig(config.content, config.command, { folder: localExamplePath(examples.folder, config.folder) });
      }

      const { command } = useDraftStore.getState().draft;
      const isCurrent = claimInputs([...new Set([...inputSettingKeys(command), ...inputSettingKeys(config.command)])]);
      const runId = await uploadRun(config.command);
      const replaces = draftUploadPaths(useDraftStore.getState().draft);
      const uploaded: UploadedConfigInput[] = [];

      const releaseUploaded = () =>
        releaseFromPendingRun(
          runId,
          uploaded.flatMap(({ stored }) => stored.map((file) => file.path)),
        );

      const loading = uploadConfigInputs(
        config,
        async (files) => (isCurrent(runId) ? uploadFiles(runId, files, replaces) : undefined),
        uploaded,
      ).then(async (complete) =>
        complete
          ? loadConfig(config.content, config.command, {
              uploaded: uploadedConfigInputs(uploaded, runId),
              isCurrent: () => isCurrent(runId),
            })
          : SUPERSEDED,
      );

      void loading.catch(releaseUploaded);

      const result = await loading;

      if (!result.loaded) {
        releaseUploaded();
      }

      return result;
    },
    [loadConfig, releaseFromPendingRun, uploadFiles, uploadRun],
  );
}

export function useDatasetLoader() {
  const navigate = useNavigate();
  const { uploadRun, uploadFiles, release, releaseFromPendingRun } = useRunUploads();

  return useCallback(
    async (examples: ExampleIndex, dataset: IndexedDataset, command: AppCommand) => {
      const inputs = indexedDatasetInputs(dataset, command);

      if (examples.origin === "local") {
        applyDraftInputs(
          inputs.map((input) => {
            const { key, value, label } = localInputAssignment(examples.folder, input);

            return { key, value, source: { label, origin: "dataset" } };
          }),
          release,
        );
        await navigate({ to: "/new" });

        return;
      }

      const isCurrent = claimInputs(inputs.map((input) => input.key));
      const files = await Promise.all(inputs.map(async (input) => fetchExampleFile(exampleDataBase(), input.path)));

      if (!isCurrent()) {
        return;
      }

      const runId = await uploadRun(command);

      if (!isCurrent(runId)) {
        return;
      }

      const replaced = new Set(inputs.map((input) => input.key));
      const stored = await uploadFiles(runId, files, draftUploadPaths(useDraftStore.getState().draft, replaced));

      if (!isCurrent(runId)) {
        releaseFromPendingRun(
          runId,
          stored.map((file) => file.path),
        );

        return;
      }

      applyDraftInputs(
        inputs.flatMap((input, index) => {
          const file = stored.slice(index, index + 1);
          const value = storedValue(file, input.list);

          return value === undefined
            ? []
            : [{ key: input.key, value, source: { ...storedSource(file), origin: "upload", run_id: runId } }];
        }),
        release,
      );
      await navigate({ to: "/new" });
    },
    [navigate, release, releaseFromPendingRun, uploadFiles, uploadRun],
  );
}

export function exampleDataBase(): string {
  return import.meta.env.BASE_URL;
}

export interface ExampleIndexQuery {
  data: ExampleIndex | undefined;
  error: Error | null;
}

async function uploadConfigInputs(
  config: IndexedConfig,
  upload: (files: readonly File[]) => Promise<StoredFile[] | undefined>,
  uploaded: UploadedConfigInput[],
): Promise<boolean> {
  for (const input of configInputs(config)) {
    const files = await Promise.all(input.paths.map(async (path) => fetchExampleFile(exampleDataBase(), path)));
    const stored = await upload(files);

    if (stored === undefined) {
      return false;
    }

    uploaded.push({ input, stored });
  }

  return true;
}

interface UploadedConfigInput {
  input: ConfigInput;
  stored: StoredFile[];
}
