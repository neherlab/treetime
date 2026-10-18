import { type JsonValue, errorMessage, type AppCommand } from "@neherlab/app-contracts";
import { useCallback, useMemo, useState } from "react";
import { useDropzone } from "react-dropzone";

import { Button } from "../ui/button";
import { Spinner } from "../ui/spinner";
import { useToastManager } from "../ui/toast";
import { dropzoneAccept } from "./fileAccept";
import { useInputActions } from "./useInputActions";

interface PathPickerProps {
  command: AppCommand;
  settingKey: string;
  title: string;
  extensions: readonly string[];
  list: boolean;
  filled: boolean;
}

export function PathPicker({ command, settingKey, title, extensions, list, filled }: PathPickerProps) {
  const actions = useInputActions(command);
  const [busy, setBusy] = useState(false);
  const toasts = useToastManager();
  const accept = useMemo(() => dropzoneAccept(extensions), [extensions]);

  const upload = useCallback(
    async (files: File[]) => {
      setBusy(true);

      try {
        await actions.addFiles(settingKey, files, list);
      } catch (error: unknown) {
        toasts.add({ title: `${title} cannot be added`, description: errorMessage(error) });
      } finally {
        setBusy(false);
      }
    },
    [actions, list, settingKey, title, toasts],
  );

  const onDropAccepted = useCallback((files: File[]) => void upload(files), [upload]);

  const dialog = useDropzone({
    multiple: list,
    noClick: true,
    noKeyboard: true,
    noDrag: true,
    ...accept,
    onDropAccepted,
  });

  const choose = useCallback(async () => {
    if (!actions.canPick) {
      dialog.open();

      return;
    }

    try {
      await actions.pick(settingKey, title, [...extensions], list);
    } catch (error: unknown) {
      toasts.add({ title: `${title} cannot be added`, description: errorMessage(error) });
    }
  }, [actions, dialog, extensions, list, settingKey, title, toasts]);

  const onChoose = useCallback(() => void choose(), [choose]);

  const onRemove = useCallback(() => actions.clear(settingKey, emptyValue(list)), [actions, list, settingKey]);

  return (
    <div className="flex shrink-0 gap-1">
      <input {...dialog.getInputProps()} aria-label={`${title} ${list ? "files" : "file"} to upload`} />
      <Button type="button" variant="outline" size="sm" onClick={onChoose} disabled={busy}>
        {busy && <Spinner />}
        {buttonText(busy, filled, list)}
      </Button>
      {filled && (
        <Button type="button" variant="ghost" size="sm" onClick={onRemove}>
          Remove
        </Button>
      )}
    </div>
  );
}

export function useFileDrop(command: AppCommand, settingKey: string, list: boolean) {
  const actions = useInputActions(command);
  const toasts = useToastManager();

  const addDropped = useCallback(
    async (files: File[]) => {
      try {
        await actions.addFiles(settingKey, files, list);
      } catch (error: unknown) {
        toasts.add({
          title: `${files.map((file) => file.name).join(", ")} cannot be added`,
          description: errorMessage(error),
        });
      }
    },
    [actions, list, settingKey, toasts],
  );

  const onDropAccepted = useCallback((files: File[]) => void addDropped(files), [addDropped]);

  return useDropzone({
    multiple: list,
    noClick: true,
    noKeyboard: true,
    noDragEventsBubbling: true,
    onDropAccepted,
  });
}

function emptyValue(list: boolean): JsonValue | undefined {
  return list ? [] : undefined;
}

function buttonText(busy: boolean, filled: boolean, list: boolean): string {
  if (busy) {
    return "Uploading...";
  }

  if (filled) {
    return "Replace";
  }

  return list ? "Choose files" : "Choose file";
}
