import type { AuspiceDocument, ResultTree } from "@neherlab/app-contracts";
import { lazy, Suspense, type ReactNode } from "react";
import TreePine from "~icons/lucide/tree-pine";

import { LoadingState } from "../components/PageShell";
import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
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
    <Suspense fallback={<LoadingState text="Loading the tree view" />}>
      <TreeWorkspace data={data} colorBy={colorBy} aside={aside} />
    </Suspense>
  );
}

export function MissingTree() {
  return (
    <Alert>
      <TreePine aria-hidden />
      <AlertTitle>No tree view</AlertTitle>
      <AlertDescription>This run wrote no Auspice tree. The output files are listed below.</AlertDescription>
    </Alert>
  );
}

export function initialColorBy(tree: ResultTree, preferred: readonly string[]): string | undefined {
  const keys = new Set(tree.colorings.map((coloring) => coloring.key));

  return preferred.find((key) => keys.has(key));
}
