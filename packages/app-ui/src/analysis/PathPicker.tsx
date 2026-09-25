import { errorMessage } from "@neherlab/app-contracts";
import type { AppCommand } from "@neherlab/app-contracts";
import { useCallback, useMemo, useRef, useState } from "react";

import type { JsonValue } from "../settings/json";
import { Button, Toast } from "../ui";
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
  const fileInput = useRef<HTMLInputElement>(null);
  const [busy, setBusy] = useState(false);
  const toasts = Toast.useToastManager();
  const accept = useMemo(() => extensions.map((extension) => `.${extension}`).join(","), [extensions]);

  const choose = useCallback(async () => {
    if (!actions.canPick) {
      fileInput.current?.click();

      return;
    }

    try {
      await actions.pick(settingKey, title, [...extensions], list);
    } catch (error: unknown) {
      toasts.add({ title: `${title} cannot be added`, description: errorMessage(error) });
    }
  }, [actions, extensions, list, settingKey, title, toasts]);

  const onChoose = useCallback(() => void choose(), [choose]);

  const onFile = useCallback(
    async (event: React.ChangeEvent<HTMLInputElement>) => {
      const file = event.target.files?.[0];
      const target = event.target;

      if (file === undefined) {
        return;
      }

      setBusy(true);

      try {
        await actions.addFile(settingKey, file, list);
      } catch (error: unknown) {
        toasts.add({ title: `${title} cannot be added`, description: errorMessage(error) });
      } finally {
        setBusy(false);
        target.value = "";
      }
    },
    [actions, list, settingKey, title, toasts],
  );

  const onFileChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => void onFile(event), [onFile]);

  const onRemove = useCallback(() => actions.clear(settingKey, emptyValue(list)), [actions, list, settingKey]);

  return (
    <div className="flex shrink-0 gap-1">
      <Button type="button" variant="outline" size="sm" onClick={onChoose} disabled={busy}>
        {buttonText(busy, filled)}
      </Button>
      {filled && (
        <Button type="button" variant="ghost" size="sm" onClick={onRemove}>
          Remove
        </Button>
      )}
      <input
        ref={fileInput}
        type="file"
        tabIndex={-1}
        aria-hidden
        accept={accept}
        onChange={onFileChange}
        className="hidden"
      />
    </div>
  );
}

export function useFileDrop(command: AppCommand, settingKey: string, list: boolean) {
  const actions = useInputActions(command);
  const [over, setOver] = useState(false);
  const toasts = Toast.useToastManager();

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

  const onDragOver = useCallback((event: React.DragEvent) => {
    event.preventDefault();
    event.stopPropagation();
    setOver(true);
  }, []);

  const onDragLeave = useCallback(() => setOver(false), []);

  const onDrop = useCallback(
    (event: React.DragEvent) => {
      event.preventDefault();
      event.stopPropagation();
      setOver(false);
      const file = event.dataTransfer.files[0];

      if (file !== undefined) {
        void addDropped(file);
      }
    },
    [addDropped],
  );

  return { over, onDragOver, onDragLeave, onDrop };
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
