import type { AppCommand, JsonValue, UiDraftSource } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";
import { useFormContext } from "react-hook-form";

import { useHost } from "../host-context";
import { baseName, type InputAssignment, pathList } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { filesForSetting, releasedPaths, type StoredFile, storedSource, storedValue } from "./fileSelection";
import { toFormValue, type FormConfig } from "./formValues";
import { useRunUploads } from "./useRunUploads";

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
  const host = useHost();
  const { getValues, setValue } = useFormContext<FormConfig>();
  const setSource = useDraftStore((state) => state.setSource);
  const { uploadRun, uploadFiles, release, releaseFromPendingRun } = useRunUploads();

  const pendingUploadPaths = useCallback(
    (key: string, runId: string): string[] => {
      const { sources, upload_run_id } = useDraftStore.getState().draft;
      const source = sources[key];
      const value = getValues(key);

      return source?.origin === "upload" && source.run_id === runId && runId === upload_run_id && value !== undefined
        ? pathList(value)
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
      const runId = await uploadRun(command);

      if (!isCurrentSelection(key, selection, runId)) {
        return;
      }

      const replaces = pendingUploadPaths(key, runId);

      const stored = await uploadFiles(runId, selected, replaces);

      if (!isCurrentSelection(key, selection, runId)) {
        releaseFromPendingRun(
          runId,
          stored.map((file) => file.path),
        );

        return;
      }

      assignStored(key, stored, list, { origin: "upload", run_id: runId });
    },
    [assignStored, command, host, pendingUploadPaths, releaseFromPendingRun, uploadFiles, uploadRun],
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
