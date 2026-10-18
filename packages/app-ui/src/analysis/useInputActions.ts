import { errorMessage, type AppCommand, type JsonValue, type UiDraftSource } from "@neherlab/app-contracts";
import { runsCreate, runsReleaseInput, runsUploadInput } from "@neherlab/app-contracts/client";
import { useCallback, useMemo } from "react";
import { useFormContext } from "react-hook-form";

import { useApiContext } from "../api/context";
import { useHost } from "../host-context";
import { baseName, type InputAssignment, pathList } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { useToastManager } from "../ui/toast";
import {
  filesForSetting,
  releasedPaths,
  type StoredFile,
  storedSource,
  storedValue,
  storeInOrder,
} from "./fileSelection";
import { fromFormValue, toFormValue, type FormConfig } from "./formValues";
import { pendingUploadRun } from "./pendingUpload";

interface InputActions {
  canPick: boolean;
  assign: (key: string, value: JsonValue, source: UiDraftSource) => void;
  assignAll: (assignments: readonly InputAssignment[], from: InputOrigin) => void;
  clear: (key: string, emptyValue: JsonValue | undefined) => void;
  addFiles: (key: string, files: readonly File[], list: boolean) => Promise<void>;
  pick: (key: string, title: string, extensions: string[], list: boolean) => Promise<void>;
}

type InputOrigin = Omit<UiDraftSource, "label" | "size">;

function isCurrentSelection(key: string, selection: symbol, runId: string): boolean {
  const { draft, selections } = useDraftStore.getState();

  return selections[key] === selection && draft.upload_run_id === runId;
}

export function useInputActions(command: AppCommand): InputActions {
  const { client } = useApiContext();
  const host = useHost();
  const { getValues, setValue } = useFormContext<FormConfig>();
  const setSource = useDraftStore((state) => state.setSource);
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

  const pendingUploadPaths = useCallback(
    (key: string, runId: string): string[] => {
      const { sources, upload_run_id } = useDraftStore.getState().draft;
      const source = sources[key];
      const value = getValues(key);

      return source?.origin === "upload" && source.run_id === runId && runId === upload_run_id && value !== undefined
        ? pathList(fromFormValue(value))
        : [];
    },
    [getValues],
  );

  const releaseReplaced = useCallback(
    (key: string, next: JsonValue | undefined) => {
      const runId = useDraftStore.getState().draft.upload_run_id;

      if (runId !== undefined) {
        release(runId, releasedPaths(pendingUploadPaths(key, runId), next));
      }
    },
    [pendingUploadPaths, release],
  );

  const assign = useCallback(
    (key: string, value: JsonValue, source: UiDraftSource) => {
      releaseReplaced(key, value);
      setValue(key, toFormValue(value), { shouldDirty: true, shouldValidate: true });
      setSource(key, source);
    },
    [releaseReplaced, setSource, setValue],
  );

  const assignAll = useCallback(
    (assignments: readonly InputAssignment[], from: InputOrigin) => {
      for (const assignment of assignments) {
        assign(assignment.key, assignment.value, { ...from, label: assignment.label });
      }
    },
    [assign],
  );

  const clear = useCallback(
    (key: string, emptyValue: JsonValue | undefined) => {
      releaseReplaced(key, emptyValue);
      setValue(key, toFormValue(emptyValue), { shouldDirty: true, shouldValidate: true });
      setSource(key, undefined);
    },
    [releaseReplaced, setSource, setValue],
  );

  const uploadRun = useCallback(async (): Promise<string> => {
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
  }, [client, command]);

  const assignStored = useCallback(
    (key: string, stored: readonly StoredFile[], list: boolean, origin: InputOrigin) => {
      const value = storedValue(stored, list);

      if (value !== undefined) {
        assign(key, value, { ...storedSource(stored), ...origin });
      }
    },
    [assign],
  );

  const addFiles = useCallback(
    async (key: string, files: readonly File[], list: boolean) => {
      const selected = filesForSetting(files, list);

      if (selected.length === 0) {
        return;
      }

      if (host !== null) {
        const stored = selected.map((file) => ({ name: file.name, path: host.pathForFile(file), size: file.size }));
        assignStored(key, stored, list, { origin: "local" });

        return;
      }

      const selection = Symbol(key);
      useDraftStore.getState().claimSelection(key, selection);
      const runId = await uploadRun();

      if (!isCurrentSelection(key, selection, runId)) {
        return;
      }

      const replaces = pendingUploadPaths(key, runId);

      const stored = await storeInOrder(
        selected,
        async (file) => {
          const { data: uploaded } = await runsUploadInput({
            client,
            path: { id: runId, name: file.name },
            query: { replaces },
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
      );

      if (!isCurrentSelection(key, selection, runId)) {
        releaseFromPendingRun(
          runId,
          stored.map((file) => file.path),
        );

        return;
      }

      assignStored(key, stored, list, { origin: "upload", run_id: runId });
    },
    [assignStored, client, host, pendingUploadPaths, releaseFromPendingRun, uploadRun],
  );

  const pick = useCallback(
    async (key: string, title: string, extensions: string[], list: boolean) => {
      if (host === null) {
        return;
      }

      const paths = await host.pickFiles({ title, extensions, multiple: list });
      const [first] = paths;

      if (first !== undefined) {
        assign(key, list ? paths : first, { label: paths.map(baseName).join(", "), origin: "local" });
      }
    },
    [assign, host],
  );

  return useMemo(
    () => ({ canPick: host !== null, assign, assignAll, clear, addFiles, pick }),
    [addFiles, assign, assignAll, clear, host, pick],
  );
}
