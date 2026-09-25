import { changeColorBy } from "auspice/src/actions/colors";
import { createStateFromQueryOrJSONs } from "auspice/src/actions/recomputeReduxState";
import { applyFilter, updateVisibleTipsAndBranchThicknesses } from "auspice/src/actions/tree";
import { CLEAN_START, SELECT_NODE } from "auspice/src/actions/types";
import { strainSymbol } from "auspice/src/util/globals";
import { useCallback, useMemo, useSyncExternalStore } from "react";

import type { AuspiceDocument } from "../results/tree";
import type { AuspiceState } from "./state";
import { createAuspiceStore, type AuspiceStore } from "./store";

export function useAuspiceStore(document: AuspiceDocument, colorBy: string | undefined): AuspiceStore {
  return useMemo(() => loadedStore(document, colorBy), [colorBy, document]);
}

export function useAuspiceSelector<T>(store: AuspiceStore, selector: (state: AuspiceState) => T): T {
  const subscribe = useCallback((listener: () => void) => store.subscribe(listener), [store]);
  const read = useCallback(() => selector(store.getState()), [selector, store]);

  return useSyncExternalStore(subscribe, read, read);
}

export function focusNode(store: AuspiceStore, name: string): void {
  const node = store.getState().tree.nodes?.find((candidate) => candidate.name === name);

  if (node === undefined) {
    return;
  }

  if (node.hasChildren) {
    store.dispatch(updateVisibleTipsAndBranchThicknesses({ root: [node.arrayIdx, undefined] }));

    return;
  }

  store.dispatch({ type: SELECT_NODE, name: node.name, idx: node.arrayIdx, isBranch: false, treeId: "LEFT" });
  store.dispatch(applyFilter("add", strainSymbol, [node.name]));
}

export function showWholeTree(store: AuspiceStore): void {
  store.dispatch(updateVisibleTipsAndBranchThicknesses({ root: [0, undefined] }));
}

function loadedStore(document: AuspiceDocument, colorBy: string | undefined): AuspiceStore {
  const store = createAuspiceStore();
  const state = createStateFromQueryOrJSONs({ json: structuredClone(document), query: {}, dispatch: store.dispatch });

  store.dispatch({ ...state, type: CLEAN_START });

  if (colorBy !== undefined) {
    store.dispatch(changeColorBy(colorBy));
  }

  return store;
}
