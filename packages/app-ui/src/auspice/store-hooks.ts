import type { AuspiceDocument } from "@neherlab/app-contracts";
import { changeColorBy } from "auspice/src/actions/colors";
import { createStateFromQueryOrJSONs } from "auspice/src/actions/recomputeReduxState";
import { applyFilter, updateVisibleTipsAndBranchThicknesses } from "auspice/src/actions/tree";
import { CLEAN_START, SELECT_NODE } from "auspice/src/actions/types";
import { genotypeSymbol, NODE_VISIBLE, strainSymbol } from "auspice/src/util/globals";
import { getTraitFromNode } from "auspice/src/util/treeMiscHelpers";
import { useCallback, useMemo, useSyncExternalStore } from "react";

import type { AuspiceControlsState, AuspiceNode, AuspiceState, AuspiceThunk } from "./state";
import { createAuspiceStore, type AuspiceStore } from "./store";

const FOCUS_MIN_CLADE_TIPS = 20;

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

  store.dispatch(showTip(node));
}

export function colorTreeBy(store: AuspiceStore, colorBy: string): void {
  store.dispatch(changeColorBy(colorBy));
}

export function showWholeTree(store: AuspiceStore): void {
  store.dispatch(updateVisibleTipsAndBranchThicknesses({ root: [0, undefined] }));
}

function showTip(tip: AuspiceNode): AuspiceThunk {
  return (dispatch, getState) => {
    const { filters } = getState().controls;

    for (const trait of Object.keys(filters)) {
      const values = activeFilterValues(filters, trait);

      if (values.length > 0 && !values.includes(String(getTraitFromNode(tip, trait)))) {
        dispatch(applyFilter("inactivate", trait, values));
      }
    }

    const { tree } = getState();
    const inView = tree.nodes?.[tree.idxOfInViewRootNode];

    if (inView === undefined || !isInClade(tip, inView)) {
      dispatch(updateVisibleTipsAndBranchThicknesses({ root: [cladeAround(tip).arrayIdx, undefined] }));
    }

    dispatch({ type: SELECT_NODE, name: tip.name, idx: tip.arrayIdx, isBranch: false, treeId: "LEFT" });
    dispatch(applyFilter("set", strainSymbol, [tip.name]));

    if (getState().tree.visibility?.[tip.arrayIdx] !== NODE_VISIBLE) {
      const genotypes = activeFilterValues(getState().controls.filters, genotypeSymbol);

      if (genotypes.length > 0) {
        dispatch(applyFilter("inactivate", genotypeSymbol, genotypes));
      }
    }
  };
}

function activeFilterValues(filters: AuspiceControlsState["filters"], trait: string | symbol): string[] {
  return (filters[trait] ?? []).filter((filter) => filter.active).map((filter) => filter.value);
}

function isInClade(node: AuspiceNode, clade: AuspiceNode): boolean {
  let current = node;

  while (current !== clade && current.parent !== current) {
    current = current.parent;
  }

  return current === clade;
}

function cladeAround(tip: AuspiceNode): AuspiceNode {
  let clade = tip.parent;

  while (clade.fullTipCount < FOCUS_MIN_CLADE_TIPS && clade.parent !== clade) {
    clade = clade.parent;
  }

  return clade;
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
