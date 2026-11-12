import { isSidebarDocked, panelSizes, renderedPanels, startPanelLayout, type ViewPanel } from "../layout";

describe("viewer start layout (auspice controlsHiddenWidth 780, twoColumnBreakpoint 1600)", () => {
  test.each([
    { width: 780, panelLayout: "full" },
    { width: 1600, panelLayout: "full" },
    { width: 1601, panelLayout: "grid" },
  ])("a viewer $width px wide starts with the $panelLayout panel layout", ({ width, panelLayout }) => {
    expect(startPanelLayout(width)).toBe(panelLayout);
  });

  test.each([
    { width: 780, docked: false },
    { width: 781, docked: true },
  ])("a viewer $width px wide docks the sidebar: $docked", ({ width, docked }) => {
    expect(isSidebarDocked(width)).toBe(docked);
  });
});

describe("panel sizes (auspice calcPanelDims: card paddings 40 and 100 px, 52 px vertical, fractions 0.88, 0.64, 0.36)", () => {
  test("one tree and the entropy chart fill the width, the tree 0.88 of the height", () => {
    expect(
      panelSizes({
        panelsWidth: 1000,
        visibleHeight: 1000,
        panelsToDisplay: ["tree", "entropy"],
        panelLayout: "full",
        twoTrees: false,
        headerStacked: false,
      }),
    ).toStrictEqual({
      tree: { width: 960, height: 802 },
      map: undefined,
      mapInGrid: false,
      entropy: { width: 960, height: 282 },
    });
  });

  test("the grid layout puts the tree and the map side by side at half width", () => {
    expect(
      panelSizes({
        panelsWidth: 1000,
        visibleHeight: 1000,
        panelsToDisplay: ["tree", "map", "entropy"],
        panelLayout: "grid",
        twoTrees: false,
        headerStacked: false,
      }),
    ).toStrictEqual({
      tree: { width: 450, height: 562 },
      map: { width: 450, height: 562 },
      mapInGrid: true,
      entropy: { width: 960, height: 282 },
    });
  });

  test("two trees take the full width with a taller title row and no map", () => {
    expect(
      panelSizes({
        panelsWidth: 1000,
        visibleHeight: 1000,
        panelsToDisplay: ["tree", "map"],
        panelLayout: "grid",
        twoTrees: true,
        headerStacked: true,
      }),
    ).toStrictEqual({ tree: { width: 960, height: 756 }, map: undefined, mapInGrid: true, entropy: undefined });
  });

  test("a stacked header row takes its height from the tree", () => {
    expect(
      panelSizes({
        panelsWidth: 1000,
        visibleHeight: 1000,
        panelsToDisplay: ["tree"],
        panelLayout: "full",
        twoTrees: false,
        headerStacked: true,
      }).tree,
    ).toStrictEqual({ width: 960, height: 776 });
  });

  test("a short view keeps the entropy chart at its smallest height of 300 px", () => {
    expect(
      panelSizes({
        panelsWidth: 1000,
        visibleHeight: 400,
        panelsToDisplay: ["tree", "entropy"],
        panelLayout: "full",
        twoTrees: false,
        headerStacked: false,
      }).entropy,
    ).toStrictEqual({ width: 960, height: 274 });
  });
});

describe("panels of the figure download (the panels of the store list that the view draws)", () => {
  test.each([
    {
      view: "a view of all panels",
      listed: ["tree", "map", "entropy"],
      sized: ["tree", "map", "entropy"],
      twoTrees: false,
      webGl: true,
      failed: [],
      panels: ["tree", "map", "entropy"],
    },
    {
      view: "a view without the entropy chart of the app",
      listed: ["tree", "map", "entropy"],
      sized: ["tree", "map"],
      twoTrees: false,
      webGl: true,
      failed: [],
      panels: ["tree", "map"],
    },
    {
      view: "a view of two trees, which hides the map",
      listed: ["tree", "map", "entropy"],
      sized: ["tree", "map", "entropy"],
      twoTrees: true,
      webGl: true,
      failed: [],
      panels: ["tree", "entropy"],
    },
    {
      view: "a view without WebGL",
      listed: ["tree", "map", "entropy"],
      sized: ["tree", "map", "entropy"],
      twoTrees: false,
      webGl: false,
      failed: [],
      panels: ["tree", "entropy"],
    },
    {
      view: "a view whose tree failed",
      listed: ["tree", "map", "entropy"],
      sized: ["tree", "map", "entropy"],
      twoTrees: false,
      webGl: true,
      failed: ["tree"],
      panels: ["map", "entropy"],
    },
    {
      view: "a dataset with panels that the view does not draw",
      listed: ["tree", "frequencies", "measurements"],
      sized: ["tree", "frequencies", "measurements"],
      twoTrees: false,
      webGl: true,
      failed: [],
      panels: ["tree"],
    },
  ] satisfies readonly ExportCase[])(
    "$view exports the panels $panels",
    ({ listed, sized, twoTrees, webGl, failed, panels }) => {
      const sizes = panelSizes({
        panelsWidth: 1000,
        visibleHeight: 1000,
        panelsToDisplay: sized,
        panelLayout: "full",
        twoTrees,
        headerStacked: false,
      });

      expect(renderedPanels({ panelsToDisplay: listed, sizes, webGl, failed: new Set(failed) })).toStrictEqual(panels);
    },
  );
});

interface ExportCase {
  view: string;
  listed: readonly string[];
  sized: readonly string[];
  twoTrees: boolean;
  webGl: boolean;
  failed: readonly ViewPanel[];
  panels: readonly string[];
}
