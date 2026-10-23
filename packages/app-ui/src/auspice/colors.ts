import { calcColorScale } from "auspice/src/util/colorScale";

import type { AuspiceMetadataState, AuspiceNode } from "./state";

const NO_CONTROLS = {};

export function traitColors(
  nodes: readonly AuspiceNode[] | null,
  metadata: AuspiceMetadataState,
  attribute: string,
): ReadonlyMap<string, string> {
  if (nodes === null || !Object.hasOwn(metadata.colorings, attribute)) {
    return new Map();
  }

  const { scale, legendValues } = calcColorScale(
    attribute,
    NO_CONTROLS,
    { nodes, visibility: null },
    undefined,
    metadata,
  );

  return new Map(legendValues.map((value) => [String(value), scale(value)]));
}
