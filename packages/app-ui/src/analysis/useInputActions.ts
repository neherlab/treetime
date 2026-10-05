import type { AppCommand, JsonValue, UiDraftOrigin } from "@neherlab/app-contracts";
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
  assign: (key: string, value: JsonValue, label: string, origin: UiDraftOrigin, size?: number) => void;
  assignAll: (assignments: readonly InputAssignment[], origin: UiDraftOrigin) => void;
  clear: (key: string, emptyValue: JsonValue | undefined) => void;
  addFile: (key: string, file: File, list: boolean) => Promise<void>;
  pick: (key: string, title: string, extensions: string[], list: boolean) => Promise<void>;
}

export function useInputActions(command: AppCommand): InputActions {
  const { client } = useApiContext();
  const host = useHost();
  const { setValue } = useFormContext<FormConfig>();
  const setSource = useDraftStore((state) => state.setSource);

  const assign = useCallback(
    (key: string, value: JsonValue, label: string, origin: UiDraftOrigin, size?: number) => {
      setValue(key, toFormValue(value), { shouldDirty: true, shouldValidate: true });
      setSource(key, size === undefined ? { label, origin } : { label, origin, size });
    },
    [setSource, setValue],
  );

  const assignAll = useCallback(
    (assignments: readonly InputAssignment[], origin: UiDraftOrigin) => {
      for (const assignment of assignments) {
        assign(assignment.key, assignment.value, assignment.label, origin);
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
      const [first] = paths;

      if (first !== undefined) {
        assign(key, list ? paths : first, paths.map(baseName).join(", "), "local");
      }
    },
    [assign, host],
  );

  return useMemo(
    () => ({ canPick: host !== null, assign, assignAll, clear, addFile, pick }),
    [addFile, assign, assignAll, clear, host, pick],
  );
}
