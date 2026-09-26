import type { AppCommand } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";
import { useFormContext } from "react-hook-form";

import { useBridge } from "../BridgeContext";
import { useLocalFiles } from "../platform";
import { baseName, type InputAssignment } from "../settings/inputs";
import type { JsonValue } from "../settings/json";
import { useDraftStore } from "../store/draft";
import type { InputOrigin } from "../store/draftSchema";
import { toFormValue, type FormConfig } from "./formValues";
import { pendingUploadRun } from "./pendingUpload";

const UPLOAD_RUN_TITLE = "Uploaded inputs";

interface InputActions {
  canPick: boolean;
  assign: (key: string, value: JsonValue, label: string, origin: InputOrigin, size: number | null) => void;
  assignAll: (assignments: readonly InputAssignment[], origin: InputOrigin) => void;
  clear: (key: string, emptyValue: JsonValue) => void;
  addFile: (key: string, file: File, list: boolean) => Promise<void>;
  pick: (key: string, title: string, extensions: string[], list: boolean) => Promise<void>;
}

export function useInputActions(command: AppCommand): InputActions {
  const bridge = useBridge();
  const localFiles = useLocalFiles();
  const { setValue } = useFormContext<FormConfig>();
  const setSource = useDraftStore((state) => state.setSource);

  const assign = useCallback(
    (key: string, value: JsonValue, label: string, origin: InputOrigin, size: number | null) => {
      setValue(key, toFormValue(value), { shouldDirty: true, shouldValidate: true });
      setSource(key, { label, origin, size });
    },
    [setSource, setValue],
  );

  const assignAll = useCallback(
    (assignments: readonly InputAssignment[], origin: InputOrigin) => {
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
    const existing = await pendingUploadRun(bridge, useDraftStore.getState().uploadRunId);

    if (existing !== undefined) {
      return existing.id;
    }

    const record = await bridge.createRun({ command, config: {}, title: UPLOAD_RUN_TITLE, defer_start: true });
    useDraftStore.getState().update({ uploadRunId: record.id });

    return record.id;
  }, [bridge, command]);

  const addFile = useCallback(
    async (key: string, file: File, list: boolean) => {
      if (localFiles !== null) {
        const path = localFiles.pathForFile(file);
        assign(key, list ? [path] : path, file.name, "local", file.size);

        return;
      }

      const uploaded = await bridge.uploadInput(await uploadRun(), file.name, file);
      assign(key, list ? [uploaded.path] : uploaded.path, file.name, "upload", uploaded.size);
    },
    [assign, bridge, localFiles, uploadRun],
  );

  const pick = useCallback(
    async (key: string, title: string, extensions: string[], list: boolean) => {
      if (localFiles === null) {
        return;
      }

      const paths = await localFiles.pickFiles({ title, extensions, multiple: list });

      if (paths.length > 0) {
        assign(key, list ? paths : (paths[0] ?? null), paths.map(baseName).join(", "), "local", null);
      }
    },
    [assign, localFiles],
  );

  return useMemo(
    () => ({ canPick: localFiles !== null, assign, assignAll, clear, addFile, pick }),
    [addFile, assign, assignAll, clear, localFiles, pick],
  );
}
