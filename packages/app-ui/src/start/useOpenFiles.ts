import { errorMessage, type UiDraftSource } from "@neherlab/app-contracts";
import { fileKinds } from "@neherlab/app-contracts/client";
import { toast } from "@neherlab/design/ui/toast";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { applyDraftInputs, type DraftInput } from "../analysis/draftInputs";
import { type StoredFile, storedSource, storedValue } from "../analysis/fileSelection";
import { useConfigLoader } from "../analysis/useConfigLoader";
import { useRunUploads } from "../analysis/useRunUploads";
import { useApiContext } from "../api/context";
import { useHost } from "../host-context";
import { commandSettings } from "../settings/catalog";
import { commandSwitchNote } from "../settings/commands";
import { folderName } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { type Placement, placeFiles } from "./filePlacement";
import { useImportRun } from "./useImportRun";

export interface OpenedFiles {
  unplaced: UnplacedName[];
}

export interface UnplacedName {
  name: string;
  reason: string;
}

export function useOpenFiles(): (files: readonly File[]) => Promise<OpenedFiles> {
  const { client } = useApiContext();
  const host = useHost();
  const navigate = useNavigate();
  const { uploadRun, uploadFiles, release } = useRunUploads();
  const loadConfig = useConfigLoader();
  const importRun = useImportRun();

  const loadYaml = useCallback(
    async (file: File) => {
      try {
        const requested = useDraftStore.getState().draft.command;
        const folder = host === null ? undefined : folderName(host.pathForFile(file));
        const result = await loadConfig(await file.text(), requested, { keepInputs: true, folder });

        toast.add(
          result.loaded
            ? { title: `Loaded ${file.name}`, description: commandSwitchNote(requested, result.command) }
            : { title: `${file.name} is not a valid config`, description: result.messages.join("; ") },
        );
      } catch (error: unknown) {
        toast.add({ title: `${file.name} cannot be loaded`, description: errorMessage(error) });
      }
    },
    [host, loadConfig],
  );

  return useCallback(
    async (files: readonly File[]): Promise<OpenedFiles> => {
      if (files.length === 0) {
        return { unplaced: [] };
      }

      const { command } = useDraftStore.getState().draft;
      const runId = host === null ? await uploadRun(command) : undefined;
      const upload = runId === undefined ? undefined : { runId, stored: await uploadFiles(runId, files, []) };

      const stored: StoredFile[] =
        upload?.stored ??
        files.map((file) => ({ name: file.name, path: host?.pathForFile(file) ?? file.name, size: file.size }));

      const { data: kinds } = await fileKinds({
        client,
        body: { files: stored.map((file) => file.path) },
        throwOnError: true,
      });

      const settings = commandSettings(command);
      const placement = placeFiles(kinds.files, settings.inputs, settings.title);
      const placed = new Set(placement.inputs.flatMap((input) => input.files));

      if (upload !== undefined) {
        release(
          upload.runId,
          upload.stored.flatMap((file, index) => (placed.has(index) ? [] : [file.path])),
        );
      }

      const assigned = assignPlacedInputs(placement, stored, runId, release);

      for (const index of placement.configs) {
        const file = files[index];

        if (file !== undefined) {
          await loadYaml(file);
        }
      }

      for (const index of placement.packages) {
        const file = files[index];

        if (file !== undefined) {
          await importRun(file);
        }
      }

      if (assigned && placement.configs.length === 0 && placement.packages.length === 0) {
        await navigate({ to: "/new" });
      }

      return {
        unplaced: placement.unplaced.flatMap(({ index, reason }) => {
          const file = files[index];

          return file === undefined ? [] : [{ name: file.name, reason }];
        }),
      };
    },
    [client, host, importRun, loadYaml, navigate, release, uploadFiles, uploadRun],
  );
}

function assignPlacedInputs(
  placement: Placement,
  stored: readonly StoredFile[],
  runId: string | undefined,
  release: (runId: string, paths: readonly string[]) => void,
): boolean {
  if (placement.inputs.length === 0) {
    return false;
  }

  if (runId !== undefined && useDraftStore.getState().draft.upload_run_id !== runId) {
    release(
      runId,
      placement.inputs.flatMap((input) => input.files.flatMap((index) => stored[index]?.path ?? [])),
    );

    return false;
  }

  const origin: Omit<UiDraftSource, "label" | "size"> =
    runId === undefined ? { origin: { kind: "local" } } : { origin: { kind: "upload" }, run_id: runId };

  const inputs: DraftInput[] = placement.inputs.flatMap((input) => {
    const chosen = input.files.flatMap((index) => stored[index] ?? []);
    const value = storedValue(chosen, input.list);

    return value === undefined ? [] : [{ key: input.key, value, source: { ...storedSource(chosen), ...origin } }];
  });

  applyDraftInputs(inputs, release);

  return true;
}
