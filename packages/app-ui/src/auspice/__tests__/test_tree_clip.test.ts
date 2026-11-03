import { setClipMask } from "auspice/src/components/tree/phyloTree/renderers";

const X_RANGE = [45, 945];

const Y_RANGE = [10, 610];

describe("tree clip area (auspice setClipMask: 5 px before the root side of the x range, 15 px after its tip side, 15 px above, 5 px below)", () => {
  test.each([
    {
      tree: "a tree without branch labels",
      branchLabels: undefined,
      orientation: [1, 1],
      area: { x: 40, y: -5, width: 920, height: 620 },
    },
    {
      tree: "a tree with branch labels in the root margin",
      branchLabels: 60,
      orientation: [1, 1],
      area: { x: -20, y: -5, width: 980, height: 620 },
    },
    {
      tree: "a mirrored tree with branch labels in the root margin",
      branchLabels: 60,
      orientation: [-1, 1],
      area: { x: 40, y: -5, width: 980, height: 620 },
    },
  ])(
    "$tree clips to $area, so that the root margin and the outermost tips stay inside",
    ({ branchLabels, orientation, area }) => {
      const tree = clipTree({ branchLabels, orientation });

      setClipMask.call(tree);

      expect(clipRect(tree)).toStrictEqual(area);
    },
  );

  test("a redraw after the branch labels change moves the clip area", () => {
    const tree = clipTree({ branchLabels: undefined, orientation: [1, 1] });
    setClipMask.call(tree);

    tree.margins.branchLabels = 60;
    setClipMask.call(tree);

    expect([tree.svg.descendants("rect").length, clipRect(tree)]).toStrictEqual([
      1,
      { x: -20, y: -5, width: 980, height: 620 },
    ]);
  });
});

function clipTree({
  branchLabels,
  orientation,
}: {
  branchLabels: number | undefined;
  orientation: readonly number[];
}): ClipTree {
  return {
    id: "main",
    xScale: { range: () => [...X_RANGE] },
    yScale: { range: () => [...Y_RANGE] },
    margins: branchLabels === undefined ? {} : { branchLabels },
    params: { orientation },
    groups: {},
    svg: new RecordingSelection("svg"),
  };
}

function clipRect(tree: ClipTree): Readonly<Record<string, string | number>> {
  const [rect] = tree.svg.descendants("rect");

  return Object.fromEntries(["x", "y", "width", "height"].map((name) => [name, rect?.attrs.get(name) ?? "missing"]));
}

interface ClipTree {
  id: string;
  xScale: { range: () => number[] };
  yScale: { range: () => number[] };
  margins: { branchLabels?: number };
  params: { orientation: readonly number[] };
  groups: { clipPath?: RecordingSelection };
  svg: RecordingSelection;
}

class RecordingSelection {
  readonly name: string;
  readonly attrs = new Map<string, string | number>();
  readonly #children: RecordingSelection[] = [];

  constructor(name: string) {
    this.name = name;
  }

  append(name: string): RecordingSelection {
    const child = new RecordingSelection(name);
    this.#children.push(child);

    return child;
  }

  select(selector: string): RecordingSelection {
    return this.descendants(selector)[0] ?? new RecordingSelection(selector);
  }

  attr(name: string, value: string | number): RecordingSelection {
    this.attrs.set(name, value);

    return this;
  }

  descendants(name: string): RecordingSelection[] {
    return this.#children.flatMap((child) => [...(child.name === name ? [child] : []), ...child.descendants(name)]);
  }
}
