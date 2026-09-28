import { createStateFromQueryOrJSONs } from "auspice/src/actions/recomputeReduxState";
import { CLEAN_START, NEW_COLORS } from "auspice/src/actions/types";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";

import { createAuspiceStore } from "../store";
import { focusNode } from "../store-hooks";

const NUC_ANNOTATION = { nuc: { start: 1, end: 3, strand: "+", type: "source" } };

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

  test.each([
    {
      annotations: "a nuc range",
      meta: { panels: ["tree", "entropy"], genome_annotations: NUC_ANNOTATION },
      panels: ["tree", "entropy"],
    },
    { annotations: "no", meta: { panels: ["tree"] }, panels: ["tree"] },
  ])(
    "a document with $annotations genome annotations displays the panels $panels (auspice entropy needs meta.genome_annotations)",
    ({ meta, panels }) => {
      const store = createAuspiceStore();

      const state = createStateFromQueryOrJSONs({
        json: auspiceDocument(meta),
        query: {},
        dispatch: store.dispatch,
      });

      store.dispatch({ ...state, type: CLEAN_START });

      expect(store.getState().controls.panelsToDisplay).toStrictEqual(panels);
    },
  );
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

function auspiceDocument(meta: Readonly<Record<string, unknown>>) {
  return {
    version: "v2",
    meta: {
      title: "TreeTime ancestral analysis",
      updated: "2026-07-19",
      colorings: [{ key: "gt", title: "Genotype", type: "categorical" }],
      display_defaults: { color_by: "gt-nuc_1" },
      ...meta,
    },
    tree: {
      name: "root",
      node_attrs: { div: 0 },
      children: [
        { name: "A", node_attrs: { div: 0.5 }, branch_attrs: { mutations: { nuc: ["A1T"] } } },
        { name: "B", node_attrs: { div: 0.1 } },
      ],
    },
    root_sequence: { nuc: "ACG" },
  };
}
