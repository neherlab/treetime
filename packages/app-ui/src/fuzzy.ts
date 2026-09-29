import uFuzzy from "@leeoniya/ufuzzy";

const MATCHER = new uFuzzy({ intraMode: 1 });

const OUT_OF_ORDER = 1;

const RANKING_THRESHOLD = 0;

export function fuzzyFilter<T>(items: readonly T[], query: string, text: (item: T) => string): T[] {
  const [indices] = MATCHER.search(
    items.map((item) => text(item)),
    query,
    OUT_OF_ORDER,
    RANKING_THRESHOLD,
  );

  if (indices === null) {
    return [...items];
  }

  const matched = new Set(indices);

  return items.filter((_, index) => matched.has(index));
}
