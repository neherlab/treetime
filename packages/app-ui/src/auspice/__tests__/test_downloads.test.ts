import { createStateFromQueryOrJSONs } from "auspice/src/actions/recomputeReduxState";
import { CLEAN_START, ENTROPY_DATA } from "auspice/src/actions/types";

import { registerDom } from "../../__tests__/dom";
import { entropyTable } from "../downloads";
import { createAuspiceStore, type AuspiceStore } from "../store";

const NUC_ANNOTATION = { nuc: { start: 1, end: 9, strand: "+", type: "source" } };

const SITE_1_ENTROPY = [{ x: 1, y: "0.693" }];

beforeAll(() => {
  registerDom();
});

describe("diversity table download (site 1 has A in one tip and T in the other: entropy ln 2, written by auspice util/entropy.js with 3 decimals)", () => {
  test("the table of an entropy chart whose bars Auspice cleared has the bars of the shown tips (auspice clears the bars of a hidden or off-screen chart, actions/entropy.js#L17-L23)", () => {
    const store = loadedStore(treeDocument());

    store.dispatch({ type: ENTROPY_DATA, data: undefined, maxYVal: 1 });

    const { entropy } = store.getState();

    expect([entropy.loaded && entropy.bars, entropyTable(store.getState())?.bars]).toStrictEqual([
      undefined,
      SITE_1_ENTROPY,
    ]);
  });

  test("the table of an entropy chart with bars has these bars", () => {
    const store = loadedStore(treeDocument());

    expect(entropyTable(store.getState())?.bars).toStrictEqual(SITE_1_ENTROPY);
  });

  test("a viewer without a dataset has no table", () => {
    expect(entropyTable(createAuspiceStore().getState())).toBeUndefined();
  });
});

function loadedStore(json: unknown): AuspiceStore {
  const store = createAuspiceStore();
  const state = createStateFromQueryOrJSONs({ json, query: {}, dispatch: store.dispatch });

  store.dispatch({ ...state, type: CLEAN_START });

  return store;
}

function treeDocument() {
  return {
    version: "v2",
    meta: {
      title: "Diversity",
      updated: "2026-11-03",
      colorings: [{ key: "country", title: "Country", type: "categorical" }],
      panels: ["tree", "entropy"],
      genome_annotations: NUC_ANNOTATION,
    },
    tree: {
      name: "root",
      node_attrs: { div: 0 },
      children: [
        { name: "A", node_attrs: { div: 0.5 }, branch_attrs: { mutations: { nuc: ["A1T"] } } },
        { name: "B", node_attrs: { div: 0.1 } },
      ],
    },
    root_sequence: { nuc: "ACGTACGTA" },
  };
}
