import { useCallback, useMemo, type ReactNode } from "react";
import { ErrorBoundary, type FallbackProps } from "react-error-boundary";

import { AuspiceTree } from "../auspice/AuspiceTree";
import type { AuspiceState } from "../auspice/state";
import type { AuspiceStore } from "../auspice/store";
import { focusNode, showWholeTree, useAuspiceSelector, useAuspiceStore } from "../auspice/store-hooks";
import type { AuspiceOutput } from "../results/load";
import type { TreeNode } from "../results/tree";

export interface TreeLink {
  focus: TreeNode | undefined;
  zoomed: boolean;
  inView: ReadonlySet<string>;
  select: (name: string) => void;
  reset: () => void;
}

export function TreeWorkspace({
  auspice,
  colorBy,
  aside,
}: {
  auspice: AuspiceOutput;
  colorBy: string | undefined;
  aside?: ((link: TreeLink) => ReactNode) | undefined;
}) {
  const store = useAuspiceStore(auspice.document, colorBy);

  return (
    <div className="grid min-w-0 gap-3.5 2xl:grid-cols-[minmax(0,1fr)_28rem]">
      <ErrorBoundary FallbackComponent={TreeFailure}>
        <AuspiceTree store={store} tips={auspice.tree.tips.length} />
      </ErrorBoundary>
      {aside !== undefined && <LinkedAside store={store} auspice={auspice} aside={aside} />}
    </div>
  );
}

function LinkedAside({
  store,
  auspice,
  aside,
}: {
  store: AuspiceStore;
  auspice: AuspiceOutput;
  aside: (link: TreeLink) => ReactNode;
}) {
  const focusName = useAuspiceSelector(store, selectFocusName);
  const zoomed = useAuspiceSelector(store, selectZoomed);
  const byName = useMemo(() => new Map(auspice.tree.nodes.map((node) => [node.name, node])), [auspice]);
  const inViewRoot = useAuspiceSelector(store, selectInViewRootName);

  const inView = useMemo(
    () => new Set(subtreeNames(byName.get(inViewRoot ?? "") ?? auspice.tree.root)),
    [auspice, byName, inViewRoot],
  );

  const select = useCallback((name: string) => focusNode(store, name), [store]);
  const reset = useCallback(() => showWholeTree(store), [store]);

  return (
    <div className="grid min-w-0 content-start gap-3.5">
      {aside({ focus: byName.get(focusName ?? ""), zoomed, inView, select, reset })}
    </div>
  );
}

function TreeFailure({ error }: FallbackProps) {
  return (
    <div role="alert" className="border-signal-danger bg-signal-danger-subtle rounded-lg border px-4 py-3.5">
      <h3 className="mb-1 font-bold">Auspice cannot draw this tree</h3>
      <p className="m-0">{error instanceof Error ? error.message : String(error)}</p>
    </div>
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

function subtreeNames(node: TreeNode): string[] {
  return [node.name, ...node.children.flatMap(subtreeNames)];
}
