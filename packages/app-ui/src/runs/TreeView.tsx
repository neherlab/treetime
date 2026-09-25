import { lazy, Suspense, type ReactNode } from "react";

import type { ResultTree } from "../results/types";
import type { JsonObject } from "../settings/json";
import type { TreeLink } from "./TreeWorkspace";

const TreeWorkspace = lazy(async () => ({ default: (await import("./TreeWorkspace")).TreeWorkspace }));

export interface TreeData {
  document: JsonObject;
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

export function initialColorBy(tree: ResultTree, preferred: readonly string[]): string | undefined {
  const keys = new Set(tree.colorings.map((coloring) => coloring.key));

  return preferred.find((key) => keys.has(key));
}
