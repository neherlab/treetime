import type { AuspiceDocument } from "@neherlab/app-contracts";
import { lazy, Suspense, type ReactNode } from "react";

import type { ResultTree } from "../results/types";
import type { TreeLink } from "./TreeWorkspace";

const TreeWorkspace = lazy(async () => ({ default: (await import("./TreeWorkspace")).TreeWorkspace }));

export interface TreeData {
  document: AuspiceDocument;
  tree: ResultTree;
}

export function TreeView({
  data,
  colorBy,
  aside,
}: {
  data: TreeData;
  colorBy: string | undefined;
  aside?: ((link: TreeLink) => ReactNode) | undefined;
}) {
  return (
    <Suspense fallback={<p className="text-ink-muted px-3.5 py-6 text-center">Loading the tree view...</p>}>
      <TreeWorkspace data={data} colorBy={colorBy} aside={aside} />
    </Suspense>
  );
}

export function MissingTree() {
  return (
    <p className="border-line bg-surface-1 text-ink-muted rounded-lg border px-4 py-3.5">
      This run wrote no Auspice tree, so the tree view is not available. The output files are listed below.
    </p>
  );
}

export function initialColorBy(tree: ResultTree, preferred: readonly string[]): string | undefined {
  const keys = new Set(tree.colorings.map((coloring) => coloring.key));

  return preferred.find((key) => keys.has(key));
}
