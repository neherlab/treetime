import { errorMessage } from "@neherlab/app-contracts";
import type { AppCommand } from "@neherlab/app-contracts";
import { useCallback, useMemo, useState } from "react";
import { type DropzoneOptions, useDropzone } from "react-dropzone";

import type { JsonValue } from "../settings/json";
import { Button } from "../ui/button";
import { Spinner } from "../ui/spinner";
import { useToastManager } from "../ui/toast";
import { useInputActions } from "./useInputActions";

const ANY_APPLICATION_TYPE = "application/*";

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
  const accept = useMemo(() => acceptExtensions(extensions), [extensions]);

  const upload = useCallback(
    async (file: File) => {
      setBusy(true);

      try {
        await actions.addFile(settingKey, file, list);
      } catch (error: unknown) {
        toasts.add({ title: `${title} cannot be added`, description: errorMessage(error) });
      } finally {
        setBusy(false);
      }
    },
    [actions, list, settingKey, title, toasts],
  );

  const onDropAccepted = useCallback(
    (files: File[]) => {
      const [file] = files;

      if (file !== undefined) {
        void upload(file);
      }
    },
    [upload],
  );

  const dialog = useDropzone({
    multiple: false,
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
      <input {...dialog.getInputProps()} />
      <Button type="button" variant="outline" size="sm" onClick={onChoose} disabled={busy}>
        {busy && <Spinner />}
        {buttonText(busy, filled)}
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
    async (file: File) => {
      try {
        await actions.addFile(settingKey, file, list);
      } catch (error: unknown) {
        toasts.add({ title: `${file.name} cannot be added`, description: errorMessage(error) });
      }
    },
    [actions, list, settingKey, toasts],
  );

  const onDropAccepted = useCallback(
    (files: File[]) => {
      const [file] = files;

      if (file !== undefined) {
        void addDropped(file);
      }
    },
    [addDropped],
  );

  return useDropzone({
    multiple: false,
    noClick: true,
    noKeyboard: true,
    noDragEventsBubbling: true,
    onDropAccepted,
  });
}

function acceptExtensions(extensions: readonly string[]): Pick<DropzoneOptions, "accept"> {
  return extensions.length === 0
    ? {}
    : { accept: { [ANY_APPLICATION_TYPE]: extensions.map((extension) => `.${extension}`) } };
}

function emptyValue(list: boolean): JsonValue {
  return list ? [] : null;
}

function buttonText(busy: boolean, filled: boolean): string {
  if (busy) {
    return "Uploading...";
  }

  return filled ? "Replace" : "Choose file";
}
