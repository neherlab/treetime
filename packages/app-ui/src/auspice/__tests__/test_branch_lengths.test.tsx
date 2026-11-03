import { within } from "@testing-library/dom";
import { createStateFromQueryOrJSONs } from "auspice/src/actions/recomputeReduxState";
import { CLEAN_START } from "auspice/src/actions/types";
import { DISTANCE } from "auspice/src/util/computeDistanceMetric";
import { createRoot, type Root } from "react-dom/client";

import { inAct, registerDom } from "../../__tests__/dom";
import { AuspiceProviders } from "../AuspiceProviders";
import { SidebarControls } from "../AuspiceSidebar";
import { figureSummary, treeFilePrefix } from "../downloads";
import { createAuspiceStore, type AuspiceStore } from "../store";
import { showSecondTree } from "../store-hooks";

const NODE_DIVS = { root: 0, A: 0.5, B: 0.1 };

const NODE_DATES = { root: 2019.9, A: 2020.04, B: 2020.2 };

const DATED_NODES = { all: ["root", "A", "B"], some: ["root", "A"], none: [] };

const NUC_AND_TWO_CDS = {
  nuc: { start: 1, end: 9, strand: "+", type: "source" },
  S: { start: 1, end: 6, strand: "+", type: "CDS" },
  E: { start: 4, end: 9, strand: "+", type: "CDS" },
};

const mounted: Root[] = [];

beforeAll(() => {
  registerDom();
});

afterEach(() => {
  for (const root of mounted.splice(0)) {
    inAct(() => root.unmount());
  }

  document.body.replaceChildren();
});

describe("auspice branch lengths", () => {
  test.each([
    { tree: "a dated tree", document: treeDocument({ dates: "all" }), measures: ["num_date", "div"] },
    { tree: "an undated tree", document: treeDocument({ dates: "none" }), measures: ["div"] },
    { tree: "a partly dated tree", document: treeDocument({ dates: "some" }), measures: ["div"] },
    {
      tree: "a tree with mutations in the genome and both CDSs",
      document: treeDocument({ dates: "none", annotations: NUC_AND_TWO_CDS, cdsMutations: { S: ["K1N"], E: ["L2F"] } }),
      measures: ["div", "mutationCount:nuc", "mutationCount:S", "mutationCount:E"],
    },
    {
      tree: "a tree without mutations in CDS E",
      document: treeDocument({ dates: "none", annotations: NUC_AND_TWO_CDS, cdsMutations: { S: ["K1N"] } }),
      measures: ["div", "mutationCount:nuc", "mutationCount:S"],
    },
  ])(
    "$tree offers the branch lengths $measures (time needs a date on every node, a mutation count needs a valid mutation on the tree)",
    ({ document, measures }) => {
      const store = loadedStore(document);

      expect(store.getState().controls.availableDistanceMeasures).toStrictEqual(measures);
    },
  );

  test("an undated tree loads by divergence without an error", () => {
    let store: AuspiceStore | undefined;

    const errors = consoleErrors(() => {
      store = loadedStore(treeDocument({ dates: "none" }));
    });

    expect([store?.getState().controls.distanceMeasure, errors]).toStrictEqual(["div", []]);
  });

  test("a second tree without dates switches both trees from time to divergence", () => {
    const store = loadedStore(treeDocument({ dates: "all" }));

    showSecondTree(store, treeDocument({ dates: "none" }), { first: "First run", second: "Second run" });

    const { controls, tree, treeToo } = store.getState();

    expect([
      controls.availableDistanceMeasures,
      controls.distanceMeasure,
      nodeDistances(tree.nodes),
      nodeDistances(treeToo.nodes),
    ]).toStrictEqual([["div"], "div", { __ROOT: 0, root: 0, A: 0.5, B: 0.1 }, { __ROOT: 0, root: 0, A: 0.5, B: 0.1 }]);
  });

  test.each([
    { tree: "a dated tree", dates: "all", summary: "Showing 2 of 2 genomes sampled between Jan 2020 and Mar 2020." },
    { tree: "a partly dated tree", dates: "some", summary: "Showing 2 of 2 genomes." },
  ] as const)("the figure summary of $tree is '$summary' (auspice datasetSummary)", ({ dates, summary }) => {
    const store = loadedStore(treeDocument({ dates }));

    expect(figureSummary(store.getState())).toStrictEqual(summary);
  });

  test.each([
    { measure: "div", prefix: "treetime-ancestral-1" },
    { measure: "num_date", prefix: "treetime-ancestral-1" },
    { measure: "mutationCount:nuc", prefix: "treetime-ancestral-1_mutations-nuc" },
    { measure: "mutationCount:S", prefix: "treetime-ancestral-1_mutations-S" },
  ])("the tree files of the branch length $measure start with $prefix", ({ measure, prefix }) => {
    expect(treeFilePrefix("treetime-ancestral-1", { distanceMeasure: measure })).toStrictEqual(prefix);
  });

  test.each([
    { tree: "a dated tree", dates: "all", controls: ["Date Range", "Branch Length"] },
    { tree: "a partly dated tree", dates: "some", controls: [] },
    { tree: "an undated tree", dates: "none", controls: [] },
  ] as const)(
    "the sidebar of $tree shows the controls $controls (one branch length leaves nothing to choose)",
    ({ dates, controls }) => {
      const store = loadedStore(treeDocument({ dates }));

      renderSidebar(store);

      const page = within(document.body);

      expect(
        [
          page.queryByText("Date Range") === null ? [] : ["Date Range"],
          page.queryByRole("combobox", { name: "Branch Length" }) === null ? [] : ["Branch Length"],
        ].flat(),
      ).toStrictEqual(controls);
    },
  );
});

function loadedStore(json: unknown): AuspiceStore {
  const store = createAuspiceStore();
  const state = createStateFromQueryOrJSONs({ json, query: {}, dispatch: store.dispatch });

  store.dispatch({ ...state, type: CLEAN_START });

  return store;
}

function renderSidebar(store: AuspiceStore): void {
  const container = document.createElement("div");
  document.body.append(container);
  const root = createRoot(container);
  mounted.push(root);

  inAct(() =>
    root.render(
      <AuspiceProviders store={store}>
        <SidebarControls store={store} entropy={false} />
      </AuspiceProviders>,
    ),
  );
}

function nodeDistances(nodes: readonly NodeWithDistance[] | null): Readonly<Record<string, unknown>> {
  return Object.fromEntries((nodes ?? []).map((node) => [node.name, node[DISTANCE]]));
}

interface NodeWithDistance {
  name: string;
  [DISTANCE]?: number;
}

function consoleErrors(run: () => void): unknown[][] {
  const original = console.error;
  const errors: unknown[][] = [];
  console.error = (...args: unknown[]) => {
    errors.push(args);
  };

  try {
    run();
  } finally {
    console.error = original;
  }

  return errors;
}

function treeDocument({
  dates,
  annotations,
  cdsMutations = {},
}: {
  dates: "all" | "some" | "none";
  annotations?: Readonly<Record<string, unknown>>;
  cdsMutations?: Readonly<Record<string, readonly string[]>>;
}): TreeDocument {
  const dated = DATED_NODES[dates];

  const document: TreeDocument = {
    version: "v2",
    meta: {
      title: "Branch lengths",
      updated: "2026-11-03",
      colorings: [{ key: "country", title: "Country", type: "categorical" }],
      panels: ["tree"],
    },
    tree: {
      name: "root",
      node_attrs: nodeAttrs("root", dated),
      children: [
        {
          name: "A",
          node_attrs: nodeAttrs("A", dated),
          branch_attrs: { mutations: { nuc: ["A1T"], ...cdsMutations } },
        },
        { name: "B", node_attrs: nodeAttrs("B", dated) },
      ],
    },
  };

  if (annotations !== undefined) {
    document.meta.genome_annotations = annotations;
    document.root_sequence = { nuc: "ACGTACGTA" };
  }

  return document;
}

function nodeAttrs(name: keyof typeof NODE_DIVS, dated: readonly string[]): NodeAttrs {
  const attrs: NodeAttrs = { div: NODE_DIVS[name] };

  if (dated.includes(name)) {
    attrs.num_date = { value: NODE_DATES[name] };
  }

  return attrs;
}

type TreeDocument = {
  version: "v2";
  meta: TreeMeta;
  tree: TreeNode;
  root_sequence?: { nuc: string };
};

interface TreeMeta {
  title: string;
  updated: string;
  colorings: readonly { key: string; title: string; type: string }[];
  panels: readonly string[];
  genome_annotations?: Readonly<Record<string, unknown>>;
}

interface TreeNode {
  name: string;
  node_attrs: NodeAttrs;
  branch_attrs?: { mutations: Readonly<Record<string, readonly string[]>> };
  children?: readonly TreeNode[];
}

interface NodeAttrs {
  div: number;
  num_date?: { value: number };
}
