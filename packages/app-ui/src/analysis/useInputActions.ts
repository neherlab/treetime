import type { AppCommand, JsonValue, UiDraftSource } from "@neherlab/app-contracts";
import { runsCreate, runsUploadInput } from "@neherlab/app-contracts/client";
import { useCallback, useMemo } from "react";
import { useFormContext } from "react-hook-form";

import { useApiContext } from "../api/context";
import { useHost } from "../host-context";
import { baseName, type InputAssignment } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { toFormValue, type FormConfig } from "./formValues";
import { pendingUploadRun } from "./pendingUpload";

interface InputActions {
  canPick: boolean;
  assign: (key: string, value: JsonValue, source: UiDraftSource) => void;
  assignAll: (assignments: readonly InputAssignment[], from: InputOrigin) => void;
  clear: (key: string, emptyValue: JsonValue | undefined) => void;
  addFile: (key: string, file: File, list: boolean) => Promise<void>;
  pick: (key: string, title: string, extensions: string[], list: boolean) => Promise<void>;
}

type InputOrigin = Omit<UiDraftSource, "label" | "size">;

export function useInputActions(command: AppCommand): InputActions {
  const { client } = useApiContext();
  const host = useHost();
  const { setValue } = useFormContext<FormConfig>();
  const setSource = useDraftStore((state) => state.setSource);

  const assign = useCallback(
    (key: string, value: JsonValue, source: UiDraftSource) => {
      setValue(key, toFormValue(value), { shouldDirty: true, shouldValidate: true });
      setSource(key, source);
    },
    [setSource, setValue],
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
      setValue(key, toFormValue(emptyValue), { shouldDirty: true, shouldValidate: true });
      setSource(key, undefined);
    },
    [setSource, setValue],
  );

  const uploadRun = useCallback(async (): Promise<string> => {
    const existing = await pendingUploadRun(client, useDraftStore.getState().draft.upload_run_id);

    if (existing !== undefined) {
      return existing.id;
    }

    const { data: record } = await runsCreate({
      client,
      body: { command, config: {}, defer_start: true },
      throwOnError: true,
    });

    useDraftStore.getState().setUploadRun(record.id);

    return record.id;
  }, [client, command]);

  const addFile = useCallback(
    async (key: string, file: File, list: boolean) => {
      if (host !== null) {
        const path = host.pathForFile(file);
        assign(key, list ? [path] : path, { label: file.name, origin: "local", size: file.size });

        return;
      }

      const runId = await uploadRun();

      const { data: uploaded } = await runsUploadInput({
        client,
        path: { id: runId, name: file.name },
        body: file,
        throwOnError: true,
      });

      assign(key, list ? [uploaded.path] : uploaded.path, {
        label: file.name,
        origin: "upload",
        size: uploaded.size,
        run_id: runId,
      });
    },
    [assign, client, host, uploadRun],
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
    () => ({ canPick: host !== null, assign, assignAll, clear, addFile, pick }),
    [addFile, assign, assignAll, clear, host, pick],
  );
}
