import { errorMessage } from "@neherlab/app-contracts";
import { toast } from "@neherlab/design/ui/toast";
import { useCallback } from "react";
import { useDropzone } from "react-dropzone";

import { type OpenedFiles, useOpenFiles } from "../start/useOpenFiles";

export function useAppFileDrop() {
  const openFiles = useOpenFiles();

  const onDropAccepted = useCallback((files: File[]) => void openDropped(openFiles, files), [openFiles]);

  return useDropzone({ noClick: true, noKeyboard: true, onDropAccepted });
}

async function openDropped(openFiles: (files: readonly File[]) => Promise<OpenedFiles>, files: readonly File[]) {
  try {
    const { unplaced } = await openFiles(files);

    if (unplaced.length > 0) {
      toast.add({
        title:
          unplaced.length === 1 ? "A dropped file was not added" : `${unplaced.length} dropped files were not added`,
        description: unplaced.map(({ name, reason }) => `${name}: ${reason}`).join("; "),
      });
    }
  } catch (error: unknown) {
    toast.add({ title: "The dropped files cannot be opened", description: errorMessage(error) });
  }
}
