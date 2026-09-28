import { CLEAN_START, NEW_COLORS } from "auspice/src/actions/types";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";

import { createAuspiceStore } from "../store";
import { focusNode } from "../store-hooks";

describe("auspice store", () => {
  beforeEach(() => {
    vi.stubGlobal("window", { innerWidth: 1280, innerHeight: 800, document: { body: { clientHeight: 800 } } });
  });

  afterEach(() => {
    vi.unstubAllGlobals();
  });

  test.each([
    { tips: 4000, skipTreeAnimation: false },
    { tips: 4001, skipTreeAnimation: true },
  ])(
    "a tree of $tips tips sets skipTreeAnimation to $skipTreeAnimation (auspice performanceFlags: above 4000 tips)",
    ({ tips, skipTreeAnimation }) => {
      const store = createAuspiceStore();

      store.dispatch(cleanStart({ tips, controls: {} }));

      expect(store.getState().controls.performanceFlags).toStrictEqual(
        new Map([["skipTreeAnimation", skipTreeAnimation]]),
      );
    },
  );

  test("a genotype color-by moves the scatterplot genotype axis to the new genotype (auspice keepScatterplotStateInSync)", () => {
    const store = createAuspiceStore();

    store.dispatch(cleanStart({ tips: 2, controls: scatterControls("scatter") }));
    store.dispatch(newGenotypeColors());

    expect(store.getState().controls.scatterVariables).toStrictEqual({
      x: "gt",
      y: "div",
      xLabel: "Genotype nuc: 2",
      xContinuous: false,
      xDomain: ["A", "C"],
      xTemporal: false,
    });
  });

  test("focusing a tip twice adds no filter", () => {
    const store = createAuspiceStore();

    store.dispatch({
      ...cleanStart({ tips: 1, controls: { filters: {} } }),
      tree: { loaded: true, nodes: [{ name: "A", hasChildren: false, arrayIdx: 0, fullTipCount: 1 }] },
    });
    focusNode(store, "A");
    focusNode(store, "A");

    expect(store.getState().controls.filters).toStrictEqual({});
  });

  test("a genotype color-by leaves the scatterplot axes unchanged in the rectangular layout", () => {
    const store = createAuspiceStore();

    store.dispatch(cleanStart({ tips: 2, controls: scatterControls("rect") }));
    store.dispatch(newGenotypeColors());

    expect(store.getState().controls.scatterVariables).toStrictEqual({ x: "gt", y: "div" });
  });
});

function cleanStart({ tips, controls }: { tips: number; controls: Readonly<Record<string, unknown>> }) {
  return {
    type: CLEAN_START,
    metadata: { loaded: true },
    tree: { loaded: true, nodes: [{ fullTipCount: tips }] },
    controls,
    entropy: {},
    measurements: {},
  };
}

function scatterControls(layout: "rect" | "scatter") {
  return {
    layout,
    colorBy: "gt-nuc_1",
    scatterVariables: { x: "gt", y: "div" },
    coloringsPresentOnTreeWithConfidence: new Set<string>(),
  };
}

function newGenotypeColors() {
  return {
    type: NEW_COLORS,
    colorBy: "gt-nuc_2",
    colorScale: { scaleType: "categorical", domain: ["A", "C"] },
  };
}
