import { errorMessage, type ResultNode, type ResultTree } from "@neherlab/app-contracts";
import { useCallback, useMemo, type ReactNode } from "react";
import { ErrorBoundary, type FallbackProps } from "react-error-boundary";

import { AuspiceTree } from "../auspice/AuspiceTree";
import type { AuspiceState } from "../auspice/state";
import type { AuspiceStore } from "../auspice/store";
import { focusNode, showWholeTree, useAuspiceSelector, useAuspiceStore } from "../auspice/store-hooks";
import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
import type { TreeData } from "./TreeView";

export interface TreeLink {
  focus: ResultNode | undefined;
  zoomed: boolean;
  inView: ReadonlySet<string>;
  select: (name: string) => void;
  reset: () => void;
}

export function TreeWorkspace({
  data,
  colorBy,
  aside,
}: {
  data: TreeData;
  colorBy: string | undefined;
  aside?: ((link: TreeLink) => ReactNode) | undefined;
}) {
  const store = useAuspiceStore(data.document, colorBy);
  const tips = useMemo(() => data.tree.nodes.filter((node) => node.children.length === 0).length, [data.tree]);

  return (
    <div className="@container min-w-0">
      <div className="grid min-w-0 grid-cols-1 gap-4 @min-[100rem]:grid-cols-[minmax(0,1fr)_28rem]">
        <ErrorBoundary FallbackComponent={TreeFailure}>
          <AuspiceTree store={store} tips={tips} />
        </ErrorBoundary>
        {aside !== undefined && <LinkedAside store={store} tree={data.tree} aside={aside} />}
      </div>
    </div>
  );
}

function LinkedAside({
  store,
  tree,
  aside,
}: {
  store: AuspiceStore;
  tree: ResultTree;
  aside: (link: TreeLink) => ReactNode;
}) {
  const focusName = useAuspiceSelector(store, selectFocusName);
  const zoomed = useAuspiceSelector(store, selectZoomed);
  const byName = useMemo(() => new Map(tree.nodes.map((node, index) => [node.name, index])), [tree]);
  const inViewRoot = useAuspiceSelector(store, selectInViewRootName);

  const inView = useMemo(
    () => new Set(subtreeNames(tree, byName.get(inViewRoot ?? "") ?? 0)),
    [tree, byName, inViewRoot],
  );

  const select = useCallback((name: string) => focusNode(store, name), [store]);
  const reset = useCallback(() => showWholeTree(store), [store]);

  return (
    <div className="grid min-w-0 content-start gap-4">
      {aside({ focus: nodeAt(tree, byName.get(focusName ?? "")), zoomed, inView, select, reset })}
    </div>
  );
}

function TreeFailure({ error }: FallbackProps) {
  return (
    <Alert variant="destructive">
      <AlertTitle>Auspice cannot draw this tree</AlertTitle>
      <AlertDescription>{errorMessage(error)}</AlertDescription>
    </Alert>
  );
}

function selectFocusName(state: AuspiceState): string | undefined {
  return state.controls.selectedNode?.name ?? selectInViewRootName(state);
}

function selectInViewRootName(state: AuspiceState): string | undefined {
  return state.tree.nodes?.[state.tree.idxOfInViewRootNode]?.name;
}

function selectZoomed(state: AuspiceState): boolean {
  return state.tree.idxOfInViewRootNode !== 0;
}

function nodeAt(tree: ResultTree, index: number | undefined): ResultNode | undefined {
  return index === undefined ? undefined : tree.nodes[index];
}

function subtreeNames(tree: ResultTree, index: number): string[] {
  const node = tree.nodes[index];

  return node === undefined ? [] : [node.name, ...node.children.flatMap((child) => subtreeNames(tree, child))];
}
