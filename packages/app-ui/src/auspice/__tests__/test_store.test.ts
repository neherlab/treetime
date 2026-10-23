import { createStateFromQueryOrJSONs } from "auspice/src/actions/recomputeReduxState";
import { applyFilter } from "auspice/src/actions/tree";
import { CLEAN_START, NEW_COLORS } from "auspice/src/actions/types";
import { NODE_VISIBLE, strainSymbol } from "auspice/src/util/globals";

import { createAuspiceStore } from "../store";
import { focusNode } from "../store-hooks";

const NUC_ANNOTATION = { nuc: { start: 1, end: 3, strand: "+", type: "source" } };

const RUNTIME_WINDOW = Object.getOwnPropertyDescriptor(globalThis, "window");

describe("auspice store", () => {
  beforeEach(() => {
    Object.defineProperty(globalThis, "window", {
      value: { innerWidth: 1280, innerHeight: 800, document: { body: { clientHeight: 800 } } },
      configurable: true,
    });
  });

  afterEach(() => {
    if (RUNTIME_WINDOW === undefined) {
      Reflect.deleteProperty(globalThis, "window");
    } else {
      Object.defineProperty(globalThis, "window", RUNTIME_WINDOW);
    }
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

  test("focusing a tip twice leaves one strain filter with that tip", () => {
    const store = loadedFocusStore();

    focusNode(store, "x1a0");
    focusNode(store, "x1a0");

    expect(store.getState().controls.filters[strainSymbol]).toStrictEqual([{ value: "x1a0", active: true }]);
  });

  test("focusing a tip outside the zoomed clade zooms to its smallest clade with at least 20 tips", () => {
    const store = loadedFocusStore();
    focusNode(store, "Y");

    focusNode(store, "x1a0");

    const { nodes, idxOfInViewRootNode } = store.getState().tree;

    expect(nodes?.[idxOfInViewRootNode]?.name).toBe("X1");
  });

  test("a tip hidden by an active trait filter becomes visible and the filter becomes inactive", () => {
    const store = loadedFocusStore();
    store.dispatch(applyFilter("add", "country", ["brazil"]));

    focusNode(store, "x1a0");

    const { tree, controls } = store.getState();
    const tip = tree.nodes?.find((node) => node.name === "x1a0");

    expect([
      tip === undefined ? undefined : tree.visibility?.[tip.arrayIdx],
      controls.filters["country"],
    ]).toStrictEqual([NODE_VISIBLE, [{ value: "brazil", active: false }]]);
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

function loadedFocusStore() {
  const store = createAuspiceStore();
  const state = createStateFromQueryOrJSONs({ json: focusDocument(), query: {}, dispatch: store.dispatch });

  store.dispatch({ ...state, type: CLEAN_START });

  return store;
}

function focusDocument() {
  return {
    version: "v2",
    meta: {
      title: "Focus",
      updated: "2026-10-24",
      colorings: [{ key: "country", title: "Country", type: "categorical" }],
      panels: ["tree"],
    },
    tree: {
      name: "root",
      node_attrs: { div: 0 },
      children: [
        clade("X", [
          clade("X1", [clade("X1a", tips("x1a", 3, "usa")), ...tips("x1", 19, "usa")]),
          ...tips("x", 8, "usa"),
        ]),
        clade("Y", tips("y", 5, "brazil")),
      ],
    },
  };
}

function clade(name: string, children: readonly unknown[]) {
  return { name, node_attrs: { div: 0.1 }, children };
}

function tips(prefix: string, count: number, country: string) {
  return Array.from({ length: count }, (_, index) => ({
    name: `${prefix}${index}`,
    node_attrs: { div: 0.2, country: { value: country } },
  }));
}
