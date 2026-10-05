import type { AppCommand, UiDraftOrigin } from "@neherlab/app-contracts";
import { runsCreate, runsUploadInput } from "@neherlab/app-contracts/client";
import { useCallback, useMemo } from "react";
import { useFormContext } from "react-hook-form";

import { useApiContext } from "../api/context";
import { useHost } from "../host-context";
import { baseName, type InputAssignment } from "../settings/inputs";
import type { JsonValue } from "../settings/json";
import { useDraftStore } from "../store/draft";
import { toFormValue, type FormConfig } from "./formValues";
import { pendingUploadRun } from "./pendingUpload";

interface InputActions {
  canPick: boolean;
  assign: (key: string, value: JsonValue, label: string, origin: UiDraftOrigin, size: number | null) => void;
  assignAll: (assignments: readonly InputAssignment[], origin: UiDraftOrigin) => void;
  clear: (key: string, emptyValue: JsonValue) => void;
  addFile: (key: string, file: File, list: boolean) => Promise<void>;
  pick: (key: string, title: string, extensions: string[], list: boolean) => Promise<void>;
}

export function useInputActions(command: AppCommand): InputActions {
  const { client } = useApiContext();
  const host = useHost();
  const { setValue } = useFormContext<FormConfig>();
  const setSource = useDraftStore((state) => state.setSource);

  const assign = useCallback(
    (key: string, value: JsonValue, label: string, origin: UiDraftOrigin, size: number | null) => {
      setValue(key, toFormValue(value), { shouldDirty: true, shouldValidate: true });
      setSource(key, { label, origin, size });
    },
    [setSource, setValue],
  );

  const assignAll = useCallback(
    (assignments: readonly InputAssignment[], origin: UiDraftOrigin) => {
      for (const assignment of assignments) {
        assign(assignment.key, assignment.value, assignment.label, origin, null);
      }
    },
    [assign],
  );

  const clear = useCallback(
    (key: string, emptyValue: JsonValue) => {
      setValue(key, toFormValue(emptyValue), { shouldDirty: true, shouldValidate: true });
      setSource(key, null);
    },
    [setSource, setValue],
  );

  const uploadRun = useCallback(async (): Promise<string> => {
    const existing = await pendingUploadRun(client, useDraftStore.getState().upload_run_id);

    if (existing !== undefined) {
      return existing.id;
    }

    const { data: record } = await runsCreate({
      client,
      body: { command, config: {}, defer_start: true },
      throwOnError: true,
    });

    useDraftStore.getState().update({ upload_run_id: record.id });

    return record.id;
  }, [client, command]);

  const addFile = useCallback(
    async (key: string, file: File, list: boolean) => {
      if (host !== null) {
        const path = host.pathForFile(file);
        assign(key, list ? [path] : path, file.name, "local", file.size);

        return;
      }

      const { data: uploaded } = await runsUploadInput({
        client,
        path: { id: await uploadRun(), name: file.name },
        body: file,
        throwOnError: true,
      });

      assign(key, list ? [uploaded.path] : uploaded.path, file.name, "upload", uploaded.size);
    },
    [assign, client, host, uploadRun],
  );

  const pick = useCallback(
    async (key: string, title: string, extensions: string[], list: boolean) => {
      if (host === null) {
        return;
      }

      const paths = await host.pickFiles({ title, extensions, multiple: list });

      if (paths.length > 0) {
        assign(key, list ? paths : (paths[0] ?? null), paths.map(baseName).join(", "), "local", null);
      }
    },
    [assign, host],
  );

  return useMemo(
    () => ({ canPick: host !== null, assign, assignAll, clear, addFile, pick }),
    [addFile, assign, assignAll, clear, host, pick],
  );
}
