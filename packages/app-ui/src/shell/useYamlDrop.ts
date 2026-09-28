import { errorMessage } from "@neherlab/app-contracts";
import { useEffect } from "react";

import { useConfigLoader } from "../analysis/useConfigLoader";
import { useDraftStore } from "../store/draft";
import { Toast } from "../ui";

const YAML_FILE = /\.ya?ml$/iu;

export function useYamlDrop() {
  const loadConfig = useConfigLoader();
  const toasts = Toast.useToastManager();

  useEffect(() => {
    function onDragOver(event: DragEvent) {
      event.preventDefault();
    }

    async function load(file: File) {
      try {
        const result = await loadConfig(await file.text(), useDraftStore.getState().command, true);

        toasts.add(
          result.loaded
            ? { title: `Loaded ${file.name}` }
            : { title: `${file.name} is not a valid config`, description: result.messages.join("; ") },
        );
      } catch (error: unknown) {
        toasts.add({
          title: `${file.name} cannot be loaded`,
          description: errorMessage(error),
        });
      }
    }

    function onDrop(event: DragEvent) {
      event.preventDefault();
      const file = event.dataTransfer?.files[0];

      if (file === undefined) {
        return;
      }

      if (YAML_FILE.test(file.name)) {
        void load(file);
      } else {
        toasts.add({ title: "Drop data files on the input slots of a new analysis" });
      }
    }

    window.addEventListener("dragover", onDragOver);
    window.addEventListener("drop", onDrop);

    return () => {
      window.removeEventListener("dragover", onDragOver);
      window.removeEventListener("drop", onDrop);
    };
  }, [loadConfig, toasts]);
}
