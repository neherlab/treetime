import { fuzzyFilter } from "../fuzzy";

const PALETTE_KINDS = ["Action", "Run", "Setting", "Example"] as const;

export function paletteItem(
  kind: PaletteKind,
  id: string,
  title: string,
  description: string,
  run: () => Promise<void> | void,
  focusId?: string,
): PaletteItem {
  return { id, kind, title, description, keywords: [kind, title, description].join(" "), run, focusId };
}

export function paletteGroups(items: readonly PaletteItem[]): PaletteGroup[] {
  const byKind = Object.groupBy(items, (item) => item.kind);

  return PALETTE_KINDS.flatMap((kind) => {
    const members = byKind[kind];

    return members === undefined ? [] : [{ kind, items: members }];
  });
}

export function matchingPaletteItems(items: readonly PaletteItem[], query: string): PaletteItem[] {
  return fuzzyFilter(items, query, (item) => item.keywords);
}

export interface PaletteItem {
  id: string;
  kind: PaletteKind;
  title: string;
  description: string;
  keywords: string;
  run: () => Promise<void> | void;
  focusId?: string | undefined;
}

export interface PaletteGroup {
  kind: PaletteKind;
  items: PaletteItem[];
}

type PaletteKind = (typeof PALETTE_KINDS)[number];
