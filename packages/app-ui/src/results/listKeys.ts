export interface Keyed<T> {
  key: string;
  item: T;
}

export function withListKeys<T>(items: readonly T[], text: (item: T) => string): Array<Keyed<T>> {
  const seen = new Map<string, number>();

  return items.map((item) => {
    const base = text(item);
    const occurrence = seen.get(base) ?? 0;
    seen.set(base, occurrence + 1);

    return { key: occurrence === 0 ? base : `${base}\u0000${occurrence}`, item };
  });
}
