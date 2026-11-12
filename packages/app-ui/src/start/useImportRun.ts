import { errorMessage, type RunRecord } from "@neherlab/app-contracts";
import { runsImport } from "@neherlab/app-contracts/client";
import { toast } from "@neherlab/design/ui/toast";
import { useNavigate } from "@tanstack/react-router";
import { useCallback } from "react";

import { useApiContext } from "../api/context";
import { useTabWorkStore } from "../notices/leaveGuard";

export function useImportRun(): (file: File) => Promise<RunRecord | undefined> {
  const { client } = useApiContext();
  const navigate = useNavigate();

  return useCallback(
    async (file: File) => {
      const endImport = useTabWorkStore.getState().startImport();

      try {
        const { data: record } = await runsImport({ client, body: file, throwOnError: true });
        await navigate({ to: "/runs/$id/results", params: { id: record.id } });

        return record;
      } catch (error: unknown) {
        toast.add({ title: `${file.name} cannot be imported`, description: errorMessage(error) });

        return undefined;
      } finally {
        endImport();
      }
    },
    [client, navigate],
  );
}
