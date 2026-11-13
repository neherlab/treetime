import type { JobId } from "@neherlab/app-contracts";
import { runsGet, runsSettings } from "@neherlab/app-contracts/client";

import { useApi } from "../api/hooks";

export function useShownRunSettings(shownRunId: JobId | undefined) {
  const enabled = shownRunId !== undefined;
  const id = shownRunId ?? "";
  const record = useApi((context) => runsGet({ ...context, path: { id } }), { enabled });

  const settings = useApi((context) => runsSettings({ ...context, path: { id } }), {
    enabled: enabled && record.data?.prepared === true,
    staleTime: Infinity,
  });

  return { record, settings };
}
