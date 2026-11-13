import { type AppCommand, errorMessage, type IndexedDataset } from "@neherlab/app-contracts";
import { examples } from "@neherlab/app-contracts/client";
import { toast } from "@neherlab/design/ui/toast";
import { useQuery } from "@tanstack/react-query";
import { useNavigate } from "@tanstack/react-router";
import { useCallback, useEffect, useMemo } from "react";

import { useApi } from "../api/hooks";
import { useHost } from "../host-context";
import { commandSettings } from "../settings/catalog";
import { normalizeConfig } from "../settings/config";
import { inputSettingKeys } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import type { QueuedExample } from "../store/draftSchema";
import { applyDraftInputs } from "./draftInputs";
import {
  type ConfigInput,
  configInputs,
  DATA_INDEX_FILE,
  indexedDatasetInputs,
  type AppConfig,
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

const EXAMPLE_COMMAND: AppCommand = "timetree";

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
    async (examples: ExampleIndex, config: AppConfig): Promise<ConfigLoadResult> => {
      if (examples.origin === "local") {
        return loadConfig(config.content, config.command, { folder: localExamplePath(examples.folder, config.folder) });
      }

      const { command } = useDraftStore.getState().draft;
      const isCurrent = claimInputs([...new Set([...inputSettingKeys(command), ...inputSettingKeys(config.command)])]);
      const runId = await uploadRun(config.command, isCurrent.selection);
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

export function useOpenExample() {
  const navigate = useNavigate();
  const { release } = useRunUploads();
  const uploadExampleInputs = useExampleInputsUpload();

  return useCallback(
    async (examples: ExampleIndex, dataset: IndexedDataset): Promise<void> => {
      const { discard } = useDraftStore.getState();
      discard("new", EXAMPLE_COMMAND);
      const { draft, load } = useDraftStore.getState();
      const specs = commandSettings(EXAMPLE_COMMAND).settings;
      load({ ...draft, config: normalizeConfig(specs, { ...draft.config, ...dataset.settings }) });
      const inputs = indexedDatasetInputs(dataset, EXAMPLE_COMMAND);

      if (examples.origin === "local") {
        applyDraftInputs(
          inputs.map((input) => {
            const { key, value, label } = localInputAssignment(examples.folder, input);

            return { key, value, source: { label, origin: { kind: "dataset" } } };
          }),
          release,
        );
        await navigate({ to: "/new" });

        return;
      }

      if (await uploadExampleInputs(inputs, EXAMPLE_COMMAND)) {
        await navigate({ to: "/new" });
      }
    },
    [navigate, release, uploadExampleInputs],
  );
}

export function useQueuedExamples(): void {
  const queued = useDraftStore((state) => state.queuedExamples);
  const command = useDraftStore((state) => state.draft.command);
  const uploadExampleInputs = useExampleInputsUpload();

  useEffect(() => {
    if (queued.length === 0) {
      return;
    }

    useDraftStore.getState().takeQueuedExamples();
    void refetchExamples(uploadExampleInputs, queued, command);
  }, [command, queued, uploadExampleInputs]);
}

async function refetchExamples(
  upload: (inputs: readonly QueuedExample[], command: AppCommand) => Promise<boolean>,
  queued: readonly QueuedExample[],
  command: AppCommand,
): Promise<void> {
  try {
    await upload(queued, command);
  } catch (error: unknown) {
    toast.add({ title: "The example files cannot be fetched again", description: errorMessage(error) });
  }
}

function useExampleInputsUpload() {
  const { uploadRun, uploadFiles, release, releaseFromPendingRun } = useRunUploads();

  return useCallback(
    async (inputs: readonly QueuedExample[], command: AppCommand): Promise<boolean> => {
      const isCurrent = claimInputs(inputs.map((input) => input.key));
      const files = await Promise.all(inputs.map(async (input) => fetchExampleFile(exampleDataBase(), input.path)));

      if (!isCurrent()) {
        return false;
      }

      const runId = await uploadRun(command, isCurrent.selection);

      if (!isCurrent(runId)) {
        return false;
      }

      const replaced = new Set(inputs.map((input) => input.key));
      const stored = await uploadFiles(runId, files, draftUploadPaths(useDraftStore.getState().draft, replaced));

      if (!isCurrent(runId)) {
        releaseFromPendingRun(
          runId,
          stored.map((file) => file.path),
        );

        return false;
      }

      applyDraftInputs(
        inputs.flatMap((input, index) => {
          const file = stored.slice(index, index + 1);
          const value = storedValue(file, input.list);

          return value === undefined
            ? []
            : [
                {
                  key: input.key,
                  value,
                  source: { ...storedSource(file), origin: { kind: "example", path: input.path }, run_id: runId },
                },
              ];
        }),
        release,
      );

      return true;
    },
    [release, releaseFromPendingRun, uploadFiles, uploadRun],
  );
}

export interface ExampleIndexQuery {
  data: ExampleIndex | undefined;
  error: Error | null;
}

async function uploadConfigInputs(
  config: AppConfig,
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

function exampleDataBase(): string {
  return import.meta.env.BASE_URL;
}
