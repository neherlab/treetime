import type { AppCommand, JsonValue, UiDraftSource } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";
import { useFormContext } from "react-hook-form";

import { useHost } from "../host-context";
import { baseName, type InputAssignment, pathList } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { applyDraftInputs } from "./draftInputs";
import {
  draftUploadPaths,
  filesForSetting,
  releasedPaths,
  type StoredFile,
  storedSource,
  storedValue,
} from "./fileSelection";
import { toFormValue, type FormConfig } from "./formValues";
import { claimInputs } from "./inputClaims";
import { useRunUploads } from "./useRunUploads";

interface InputActions {
  canPick: boolean;
  assign: (key: string, value: JsonValue, source: UiDraftSource) => void;
  assignAll: (assignments: readonly InputAssignment[], from: InputOrigin) => void;
  clear: (key: string, emptyValue: JsonValue | undefined) => void;
  addFiles: (key: string, files: readonly File[], list: boolean) => Promise<void>;
  addPendingFiles: (key: string, pending: Promise<readonly File[]>, list: boolean) => Promise<void>;
  pick: (key: string, title: string, extensions: string[], list: boolean) => Promise<void>;
}

type InputOrigin = Omit<UiDraftSource, "label" | "size">;

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

  const assignPending = useCallback(
    (key: string, stored: readonly StoredFile[], list: boolean, origin: InputOrigin, epoch: number) => {
      const value = storedValue(stored, list);

      if (value === undefined) {
        return;
      }

      const source = { ...storedSource(stored), ...origin };

      if (useDraftStore.getState().epoch === epoch) {
        assign(key, value, source);
      } else {
        applyDraftInputs([{ key, value, source }], release);
      }
    },
    [assign, release],
  );

  const addPendingFiles = useCallback(
    async (key: string, pending: Promise<readonly File[]>, list: boolean) => {
      const { epoch } = useDraftStore.getState();
      const isCurrent = claimInputs([key]);
      const selected = filesForSetting(await pending, list);

      if (host !== null) {
        const stored = selected.map((file) => ({ name: file.name, path: host.pathForFile(file), size: file.size }));

        if (isCurrent()) {
          assignPending(key, stored, list, { origin: "local" }, epoch);
        }

        return;
      }

      if (selected.length === 0 || !isCurrent()) {
        return;
      }

      const runId = await uploadRun(command);

      if (!isCurrent(runId)) {
        return;
      }

      const { draft, epoch: current } = useDraftStore.getState();
      const replaces = current === epoch ? pendingUploadPaths(key, runId) : draftUploadPaths(draft, new Set([key]));
      const stored = await uploadFiles(runId, selected, replaces);

      if (!isCurrent(runId)) {
        releaseFromPendingRun(
          runId,
          stored.map((file) => file.path),
        );

        return;
      }

      assignPending(key, stored, list, { origin: "upload", run_id: runId }, epoch);
    },
    [assignPending, command, host, pendingUploadPaths, releaseFromPendingRun, uploadFiles, uploadRun],
  );

  const addFiles = useCallback(
    async (key: string, files: readonly File[], list: boolean) => {
      if (filesForSetting(files, list).length > 0) {
        await addPendingFiles(key, Promise.resolve(files), list);
      }
    },
    [addPendingFiles],
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
    () => ({ canPick: host !== null, assign, assignAll, clear, addFiles, addPendingFiles, pick }),
    [addFiles, addPendingFiles, assign, assignAll, clear, host, pick],
  );
}
