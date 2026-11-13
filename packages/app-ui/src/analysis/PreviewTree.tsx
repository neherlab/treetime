import type { AuspiceDocument } from "@neherlab/app-contracts";
import { version as getVersion } from "@neherlab/app-contracts/client";
import { useMemo } from "react";
import { ErrorBoundary } from "react-error-boundary";

import { useApi } from "../api/hooks";
import { AuspiceTree } from "../auspice/AuspiceTree";
import type { FigureRun } from "../auspice/downloads";
import { useAuspiceStore } from "../auspice/store-hooks";
import { TreeFailure } from "../runs/TreeWorkspace";

export function PreviewTree({ document, run }: { document: AuspiceDocument; run: FigureRun }) {
  const store = useAuspiceStore(document, undefined);
  const runs = useMemo(() => [run] as const, [run]);
  const { data: version } = useApi((context) => getVersion(context), { staleTime: Infinity });

  return (
    <ErrorBoundary FallbackComponent={TreeFailure}>
      <AuspiceTree store={store} runs={runs} version={version?.version} entropy={false} />
    </ErrorBoundary>
  );
}
