const DECIMAL = /^[+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/iu;

export function parseNumber(text: string): number | string {
  const trimmed = text.trim();
  const number = DECIMAL.test(trimmed) ? Number(trimmed) : Number.NaN;

  return Number.isFinite(number) ? number : text;
}
