import type { ResultColoring } from "./types";

const MUTED_PALETTE = [
  "#332288",
  "#88ccee",
  "#44aa99",
  "#117733",
  "#999933",
  "#ddcc77",
  "#cc6677",
  "#882255",
  "#aa4499",
] as const;

const TREETIME_CATEGORICAL = new Set(["gt", "bad_branch"]);

export type ColorScales = ReadonlyMap<string, ReadonlyArray<readonly [string, string]>>;

export function mutedColorScales(colorings: readonly ResultColoring[]): ColorScales {
  return new Map(
    colorings.flatMap((coloring) =>
      coloring.kind !== "categorical" ||
      TREETIME_CATEGORICAL.has(coloring.key) ||
      coloring.states.length === 0 ||
      coloring.states.length > MUTED_PALETTE.length
        ? []
        : [
            [
              coloring.key,
              coloring.states.map((state, index) => [state, MUTED_PALETTE[index] ?? MUTED_PALETTE[0]] as const),
            ] as const,
          ],
    ),
  );
}
