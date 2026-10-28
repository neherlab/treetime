import { type JsonValue, errorMessage, type AppCommand, type IndexedFile } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { Spinner } from "@neherlab/design/ui/spinner";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useCallback, useMemo, useState } from "react";
import { useDropzone } from "react-dropzone";

import { type ExampleIndex, fetchExampleFile, localExamplePath } from "./exampleData";
import { ExampleFilesDialog } from "./ExampleFilesDialog";
import { dropzoneAccept } from "./fileAccept";
import { exampleDataBase, useExampleIndex } from "./useExampleData";
import { useInputActions } from "./useInputActions";

const NO_FILES: readonly IndexedFile[] = [];

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
  const { data: examples } = useExampleIndex();
  const [examplesOpen, setExamplesOpen] = useState(false);

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

  const addExampleFile = useCallback(
    async (example: ExampleIndex, file: IndexedFile) => {
      if (example.origin === "local") {
        const path = localExamplePath(example.folder, file.path);
        actions.assign(settingKey, list ? [path] : path, { label: file.path, origin: "dataset", size: file.size });

        return;
      }

      setBusy(true);

      try {
        await actions.addPendingFiles(
          settingKey,
          fetchExampleFile(exampleDataBase(), file.path).then((fetched) => [fetched]),
          list,
        );
      } catch (error: unknown) {
        toasts.add({ title: `${file.path} cannot be added`, description: errorMessage(error) });
      } finally {
        setBusy(false);
      }
    },
    [actions, list, settingKey, toasts],
  );

  const chooseExample = useCallback(
    (file: IndexedFile) => {
      setExamplesOpen(false);

      if (examples !== undefined) {
        void addExampleFile(examples, file);
      }
    },
    [examples, addExampleFile],
  );

  const onOpenExamples = useCallback(() => setExamplesOpen(true), []);
  const exampleFiles = examples?.index?.files ?? NO_FILES;

  return (
    <div className="flex shrink-0 gap-1">
      <input {...dialog.getInputProps()} aria-label={`${title} ${list ? "files" : "file"} to upload`} />
      <Button type="button" variant="outline" size="sm" onClick={onChoose} disabled={busy}>
        {busy && <Spinner />}
        {buttonText(busy, filled, list)}
      </Button>
      {exampleFiles.length > 0 && (
        <>
          <Button type="button" variant="ghost" size="sm" onClick={onOpenExamples} disabled={busy}>
            Example files
          </Button>
          <ExampleFilesDialog
            open={examplesOpen}
            onOpenChange={setExamplesOpen}
            title={title}
            files={exampleFiles}
            choose={chooseExample}
          />
        </>
      )}
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
