import type { AuspiceMetadataState, AuspiceNode } from "../../state";

export declare function calcColorScale(
  colorBy: string,
  controls: { colorScale?: { version: number } },
  tree: { nodes: readonly AuspiceNode[]; visibility: readonly number[] | null },
  treeToo: undefined,
  metadata: AuspiceMetadataState,
): AuspiceColorScale;

interface AuspiceColorScale {
  scale: (value: string | number) => string;
  legendValues: readonly (string | number)[];
}
