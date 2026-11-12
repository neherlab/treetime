import { fuzzyFilter } from "../fuzzy";
import type { AuspiceColoring } from "./state";

export const GENOTYPE_COLORING = "gt";

const TIP_MATCH_LIMIT = 50;

const GENOTYPE_COLOR_BY = /^gt-(.+)_(\d+(?:,\d+)*)$/u;

const TREE_LAYOUTS = [
  { value: "rect", label: "Rectangular" },
  { value: "radial", label: "Radial" },
  { value: "unrooted", label: "Unrooted" },
] as const;

const CLOCK_LAYOUT = { value: "clock", label: "Clock" } as const;

export interface ToolbarOption {
  value: string;
  label: string;
}

export function colorOptions(colorings: Readonly<Record<string, AuspiceColoring>>, colorBy: string): ToolbarOption[] {
  const options = Object.entries(colorings).flatMap(([key, coloring]) =>
    key === GENOTYPE_COLORING ? [] : [{ value: key, label: coloring.title ?? key }],
  );

  return GENOTYPE_COLOR_BY.test(colorBy) && !options.some((option) => option.value === colorBy)
    ? [{ value: colorBy, label: genotypeLabel(colorBy) }, ...options]
    : options;
}

export function layoutOptions(distanceMeasures: readonly string[]): ToolbarOption[] {
  const clock = distanceMeasures.includes("div") && distanceMeasures.includes("num_date");

  return clock ? [...TREE_LAYOUTS, CLOCK_LAYOUT] : [...TREE_LAYOUTS];
}

export function tipNames(nodes: readonly { name: string; hasChildren: boolean }[] | null): string[] {
  return (nodes ?? []).flatMap((node) => (node.hasChildren ? [] : [node.name])).toSorted();
}

export function tipMatches(tips: readonly string[], query: string): string[] {
  const trimmed = query.trim();

  return trimmed === "" ? [] : fuzzyFilter(tips, trimmed, (tip) => tip).slice(0, TIP_MATCH_LIMIT);
}

export function genotypeColorBy(position: number): string {
  return `${GENOTYPE_COLORING}-nuc_${position}`;
}

function genotypeLabel(colorBy: string): string {
  const match = GENOTYPE_COLOR_BY.exec(colorBy);
  const gene = match?.[1];
  const positions = match?.[2];

  return gene === undefined || positions === undefined
    ? colorBy
    : `Genotype ${gene} ${positions.replaceAll(",", ", ")}`;
}
