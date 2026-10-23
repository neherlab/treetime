import { changeColorBy } from "auspice/src/actions/colors";
import { createStateFromQueryOrJSONs } from "auspice/src/actions/recomputeReduxState";
import { CLEAN_START } from "auspice/src/actions/types";

import { traitColors } from "../colors";
import { createAuspiceStore } from "../store";

const RUNTIME_WINDOW = Object.getOwnPropertyDescriptor(globalThis, "window");

describe("trait colors", () => {
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
    { name: "a trait with more states than the app palette", states: 12, scale: undefined },
    {
      name: "a trait with the colors of the file",
      states: 3,
      scale: [
        ["s0", "#332288"],
        ["s1", "#88ccee"],
        ["s2", "#44aa99"],
      ],
    },
  ])("$name gives each state the color of its tips in the tree", ({ states, scale }) => {
    const store = createAuspiceStore();

    const state = createStateFromQueryOrJSONs({
      json: traitDocument(states, scale),
      query: {},
      dispatch: store.dispatch,
    });

    store.dispatch({ ...state, type: CLEAN_START });
    store.dispatch(changeColorBy("country"));

    const { tree, metadata } = store.getState();

    const treeColors = new Map(
      (tree.nodes ?? [])
        .filter((node) => !node.hasChildren)
        .map((node) => [`s${node.name.slice(3)}`, tree.nodeColors[node.arrayIdx]]),
    );

    expect([traitColors(tree.nodes, metadata, "country"), new Set(treeColors.values()).size]).toStrictEqual([
      treeColors,
      states,
    ]);
  });

  test("a coloring the file does not have gives no colors", () => {
    const store = createAuspiceStore();

    const state = createStateFromQueryOrJSONs({
      json: traitDocument(2, undefined),
      query: {},
      dispatch: store.dispatch,
    });

    store.dispatch({ ...state, type: CLEAN_START });

    const { tree, metadata } = store.getState();

    expect(traitColors(tree.nodes, metadata, "region")).toStrictEqual(new Map());
  });
});

function traitDocument(states: number, scale: readonly (readonly string[])[] | undefined) {
  const coloring = { key: "country", title: "Country", type: "categorical" };

  return {
    version: "v2",
    meta: {
      title: "Trait colors",
      updated: "2026-10-24",
      colorings: [scale === undefined ? coloring : { ...coloring, scale }],
      panels: ["tree"],
    },
    tree: {
      name: "root",
      node_attrs: { div: 0 },
      children: Array.from({ length: states }, (_, index) => ({
        name: `tip${index}`,
        node_attrs: { div: 0.1, country: { value: `s${index}` } },
      })),
    },
  };
}
