import type { AppCommand, IndexedConfig, IndexedDataset } from "@neherlab/app-contracts";
import { examples } from "@neherlab/app-contracts/client";
import { useQuery } from "@tanstack/react-query";
import { useCallback, useMemo } from "react";

import { useApi } from "../api/hooks";
import { useHost } from "../host-context";
import { useDraftStore } from "../store/draft";
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
import type { StoredFile } from "./fileSelection";
import { type ConfigLoadResult, useConfigLoader } from "./useConfigLoader";
import { useInputActions } from "./useInputActions";
import { useRunUploads } from "./useRunUploads";

const DISCARDED: ConfigLoadResult = {
  loaded: false,
  messages: ["The form was reset while the example files were uploaded"],
};

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

      const runId = await uploadRun(config.command);
      const uploaded: UploadedConfigInput[] = [];

      const releaseUploaded = () =>
        releaseFromPendingRun(
          runId,
          uploaded.flatMap(({ stored }) => stored.map((file) => file.path)),
        );

      const loading = uploadConfigInputs(config, async (files) => uploadFiles(runId, files, []), uploaded).then(
        async () =>
          useDraftStore.getState().draft.upload_run_id === runId
            ? loadConfig(config.content, config.command, { uploaded: uploadedConfigInputs(uploaded, runId) })
            : DISCARDED,
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

export function useDatasetLoader(command: AppCommand) {
  const { addFiles, assignAll } = useInputActions(command);

  return useCallback(
    async (examples: ExampleIndex, dataset: IndexedDataset) => {
      const inputs = indexedDatasetInputs(dataset, command);

      if (examples.origin === "local") {
        assignAll(
          inputs.map((input) => localInputAssignment(examples.folder, input)),
          { origin: "dataset" },
        );

        return;
      }

      const fetched = await Promise.all(
        inputs.map(async (input) => ({ input, file: await fetchExampleFile(exampleDataBase(), input.path) })),
      );

      await Promise.all(fetched.map(async ({ input, file }) => addFiles(input.key, [file], input.list)));
    },
    [addFiles, assignAll, command],
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
  upload: (files: readonly File[]) => Promise<StoredFile[]>,
  uploaded: UploadedConfigInput[],
): Promise<void> {
  for (const input of configInputs(config)) {
    const files = await Promise.all(input.paths.map(async (path) => fetchExampleFile(exampleDataBase(), path)));

    uploaded.push({ input, stored: await upload(files) });
  }
}

interface UploadedConfigInput {
  input: ConfigInput;
  stored: StoredFile[];
}
