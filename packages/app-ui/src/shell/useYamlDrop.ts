import { errorMessage } from "@neherlab/app-contracts";
import { useCallback } from "react";
import { useDropzone, type FileRejection } from "react-dropzone";

import { useConfigLoader } from "../analysis/useConfigLoader";
import { useHost } from "../host-context";
import { commandSwitchNote } from "../settings/commands";
import { folderName } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { useToastManager } from "../ui/toast";

const YAML_ACCEPT = { "application/yaml": [".yaml", ".yml"] };

export function useYamlDrop() {
  const loadConfig = useConfigLoader();
  const host = useHost();
  const toasts = useToastManager();

  const load = useCallback(
    async (file: File) => {
      try {
        const requested = useDraftStore.getState().command;
        const folder = host === null ? undefined : folderName(host.pathForFile(file));
        const result = await loadConfig(await file.text(), requested, true, folder);

        toasts.add(
          result.loaded
            ? { title: `Loaded ${file.name}`, description: commandSwitchNote(requested, result.command) }
            : { title: `${file.name} is not a valid config`, description: result.messages.join("; ") },
        );
      } catch (error: unknown) {
        toasts.add({ title: `${file.name} cannot be loaded`, description: errorMessage(error) });
      }
    },
    [host, loadConfig, toasts],
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
