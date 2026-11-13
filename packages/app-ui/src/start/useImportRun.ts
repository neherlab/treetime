import { errorMessage, type RunRecord } from "@neherlab/app-contracts";
import { runsImport, runsImportFile } from "@neherlab/app-contracts/client";
import { toast } from "@neherlab/design/ui/toast";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import { useHost } from "../host-context";
import { useTabWorkStore } from "../notices/leaveGuard";

export function useImportRun(): (file: File) => Promise<RunRecord | undefined> {
  const { client } = useApiContext();
  const navigate = useNavigate();
  const host = useHost();

  return useCallback(
    async (file: File) => {
      const endImport = useTabWorkStore.getState().startImport();

      try {
        const { data: record } =
          host === null
            ? await runsImport({ client, body: file, throwOnError: true })
            : await runsImportFile({ client, body: { path: host.pathForFile(file) }, throwOnError: true });

        await navigate({ to: "/runs/$id/results", params: { id: record.id } });

        return record;
      } catch (error: unknown) {
        toast.add({ title: `${file.name} cannot be imported`, description: errorMessage(error) });

        return undefined;
      } finally {
        endImport();
      }
    },
    [client, host, navigate],
  );
}
