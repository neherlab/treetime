export function wordMatcher(query: string): (text: string) => boolean {
  const words = query
    .toLowerCase()
    .split(/\s+/u)
    .filter((word) => word !== "");

  return (text) => {
    const haystack = text.toLowerCase();

    return words.every((word) => haystack.includes(word));
  };
}
