import { isSidebarDocked, panelSizes, startLayout } from "../layout";

describe("viewer start layout (auspice controlsHiddenWidth 780, twoColumnBreakpoint 1600)", () => {
  test.each([
    { width: 780, sidebarOpen: false, panelLayout: "full" },
    { width: 781, sidebarOpen: true, panelLayout: "full" },
    { width: 1600, sidebarOpen: true, panelLayout: "full" },
    { width: 1601, sidebarOpen: true, panelLayout: "grid" },
  ])(
    "a viewer $width px wide starts with the sidebar open: $sidebarOpen, layout $panelLayout",
    ({ width, ...expected }) => {
      expect(startLayout(width)).toStrictEqual(expected);
    },
  );

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
