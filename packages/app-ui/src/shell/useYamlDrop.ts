import { errorMessage } from "@neherlab/app-contracts";
import { useCallback } from "react";
import { useDropzone, type FileRejection } from "react-dropzone";

import { useConfigLoader } from "../analysis/useConfigLoader";
import { commandSwitchNote } from "../settings/commands";
import { useDraftStore } from "../store/draft";
import { useToastManager } from "../ui/toast";

const YAML_ACCEPT = { "application/yaml": [".yaml", ".yml"] };

export function useYamlDrop() {
  const loadConfig = useConfigLoader();
  const toasts = useToastManager();

  const load = useCallback(
    async (file: File) => {
      try {
        const requested = useDraftStore.getState().command;
        const result = await loadConfig(await file.text(), requested, true);

        toasts.add(
          result.loaded
            ? { title: `Loaded ${file.name}`, description: commandSwitchNote(requested, result.command) }
            : { title: `${file.name} is not a valid config`, description: result.messages.join("; ") },
        );
      } catch (error: unknown) {
        toasts.add({ title: `${file.name} cannot be loaded`, description: errorMessage(error) });
      }
    },
    [loadConfig, toasts],
  );

  const onDropAccepted = useCallback(
    (files: File[]) => {
      const [file] = files;

      if (file !== undefined) {
        void load(file);
      }
    },
    [load],
  );

  const onDropRejected = useCallback(
    (rejections: FileRejection[]) => {
      if (rejections.length > 0) {
        toasts.add({ title: "Drop data files on the input slots of a new analysis" });
      }
    },
    [toasts],
  );

  return useDropzone({
    accept: YAML_ACCEPT,
    multiple: false,
    noClick: true,
    noKeyboard: true,
    onDropAccepted,
    onDropRejected,
  });
}
