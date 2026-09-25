import { lazy, Suspense, type ReactNode } from "react";

import type { AuspiceOutput } from "../results/load";
import type { TreeLink } from "./TreeWorkspace";

const TreeWorkspace = lazy(async () => ({ default: (await import("./TreeWorkspace")).TreeWorkspace }));

export function TreeView({
  auspice,
  colorBy,
  aside,
}: {
  auspice: AuspiceOutput;
  colorBy: string | undefined;
  aside?: ((link: TreeLink) => ReactNode) | undefined;
}) {
  return (
    <Suspense fallback={<p className="text-ink-muted px-3.5 py-6 text-center">Loading the tree view...</p>}>
      <TreeWorkspace auspice={auspice} colorBy={colorBy} aside={aside} />
    </Suspense>
  );
}
