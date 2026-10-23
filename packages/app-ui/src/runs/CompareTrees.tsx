import type { AuspiceDocument } from "@neherlab/app-contracts";
import { version as getVersion } from "@neherlab/app-contracts/client";
import { useMemo } from "react";
import { ErrorBoundary } from "react-error-boundary";

import { useApi } from "../api/hooks";
import { AuspiceTree } from "../auspice/AuspiceTree";
import type { FigureRun } from "../auspice/downloads";
import { useTwoTreeStore } from "../auspice/store-hooks";
import { TreeFailure } from "./TreeWorkspace";

export function CompareTrees({ first, second }: { first: RunTree; second: RunTree }) {
  const names = useMemo(() => ({ first: first.run.title, second: second.run.title }), [first.run, second.run]);
  const store = useTwoTreeStore(first.document, second.document, names);
  const runs = useMemo(() => [first.run, second.run] as const, [first.run, second.run]);
  const { data: version } = useApi((context) => getVersion(context), { staleTime: Infinity });

  return (
    <ErrorBoundary FallbackComponent={TreeFailure}>
      <AuspiceTree key={`${first.run.id}:${second.run.id}`} store={store} runs={runs} version={version?.version} />
    </ErrorBoundary>
  );
}

export interface RunTree {
  run: FigureRun;
  document: AuspiceDocument;
}
