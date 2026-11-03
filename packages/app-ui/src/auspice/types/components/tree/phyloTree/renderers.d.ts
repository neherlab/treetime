export declare function setClipMask(this: ClipMaskTree): void;

interface ClipMaskTree {
  id: string | number;
  xScale: { range: () => number[] };
  yScale: { range: () => number[] };
  margins: { branchLabels?: number };
  params: { orientation: readonly number[] };
  groups: { clipPath?: ClipMaskSelection };
  svg: ClipMaskSelection;
}

interface ClipMaskSelection {
  append: (name: string) => ClipMaskSelection;
  select: (selector: string) => ClipMaskSelection;
  attr: (name: string, value: string | number) => ClipMaskSelection;
}
