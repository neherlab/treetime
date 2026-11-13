import { type DetectedFile, errorMessage, type UiDraftSource } from "@neherlab/app-contracts";
import { type ApiClient, fileKinds, fileRoles } from "@neherlab/app-contracts/client";
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
import { detectedFiles, fileHead, type RoledFile, withRoles } from "./openedFiles";
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

      const { data: roles } = await fileRoles({
        client,
        body: { files: await Promise.all(files.map(fileHead)) },
        throwOnError: true,
      });

      const opened = withRoles(files, roles.files);
      const inputs = opened.filter(({ role }) => role === "input");
      const { command } = useDraftStore.getState().draft;
      const runId = host === null && inputs.length > 0 ? await uploadRun(command) : undefined;

      const read = await readInputs(client, inputs, (files) =>
        runId === undefined
          ? files.map((file) => ({ name: file.name, path: host?.pathForFile(file) ?? file.name, size: file.size }))
          : uploadFiles(runId, files, []),
      );

      const settings = commandSettings(command);
      const placement = placeFiles(detectedFiles(opened, read.kinds), settings.inputs, settings.title);
      const placed = new Set(placement.inputs.flatMap((input) => input.files));

      if (runId !== undefined) {
        release(
          runId,
          inputs.flatMap(({ index }) => (placed.has(index) ? [] : (read.stored.get(index)?.path ?? []))),
        );
      }

      const assigned = assignPlacedInputs(placement, read.stored, runId, release);

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

interface ReadInputs {
  stored: ReadonlyMap<number, StoredFile | undefined>;
  kinds: DetectedFile[];
}

async function readInputs(
  client: ApiClient,
  inputs: ReadonlyArray<RoledFile<File>>,
  store: (files: readonly File[]) => StoredFile[] | Promise<StoredFile[]>,
): Promise<ReadInputs> {
  if (inputs.length === 0) {
    return { stored: new Map(), kinds: [] };
  }

  const stored = await store(inputs.map(({ file }) => file));

  const { data: kinds } = await fileKinds({
    client,
    body: { files: stored.map((file) => file.path) },
    throwOnError: true,
  });

  return { stored: new Map(inputs.map(({ index }, position) => [index, stored[position]])), kinds: kinds.files };
}

function assignPlacedInputs(
  placement: Placement,
  stored: ReadonlyMap<number, StoredFile | undefined>,
  runId: string | undefined,
  release: (runId: string, paths: readonly string[]) => void,
): boolean {
  if (placement.inputs.length === 0) {
    return false;
  }

  if (runId !== undefined && useDraftStore.getState().draft.upload_run_id !== runId) {
    release(
      runId,
      placement.inputs.flatMap((input) => input.files.flatMap((index) => stored.get(index)?.path ?? [])),
    );

    return false;
  }

  const origin: Omit<UiDraftSource, "label" | "size"> =
    runId === undefined ? { origin: { kind: "local" } } : { origin: { kind: "upload" }, run_id: runId };

  const inputs: DraftInput[] = placement.inputs.flatMap((input) => {
    const chosen = input.files.flatMap((index) => stored.get(index) ?? []);
    const value = storedValue(chosen, input.list);

    return value === undefined ? [] : [{ key: input.key, value, source: { ...storedSource(chosen), ...origin } }];
  });

  applyDraftInputs(inputs, release);

  return true;
}
