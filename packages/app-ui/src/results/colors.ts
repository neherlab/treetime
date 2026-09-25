import { parseAuspiceJson, type AuspiceDocument, type ResultTree } from "./tree";

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

function categoricalScale(states: readonly string[]): Array<[string, string]> | undefined {
  const distinct = [...new Set(states)].toSorted((left, right) => left.localeCompare(right));

  return distinct.length <= MUTED_PALETTE.length
    ? distinct.map((state, index) => [state, MUTED_PALETTE[index] ?? MUTED_PALETTE[0]])
    : undefined;
}

export function mutedAuspiceDocument(text: string, tree: ResultTree): AuspiceDocument {
  const document = parseAuspiceJson(text);

  const scales = new Map(
    (document.meta.colorings ?? []).flatMap((coloring) => {
      if (coloring.type !== "categorical" || TREETIME_CATEGORICAL.has(coloring.key)) {
        return [];
      }

      const scale = categoricalScale(tree.nodes.flatMap((node) => node.traits.get(coloring.key)?.value ?? []));

      return scale === undefined ? [] : [[coloring, scale] as const];
    }),
  );

  for (const [coloring, scale] of scales) {
    coloring["scale"] = scale;
  }

  return document;
}
