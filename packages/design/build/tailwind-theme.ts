import { ALIAS_COLORS, BASE_COLORS, CORNERS, TYPEFACES } from "../src/tokens/tokens";

const TAILWIND_RADII = {
  sm: "sm",
  md: "md",
  lg: "lg",
  xl: "xl",
  "2xl": "xl",
  "3xl": "xl",
  "4xl": "xl",
} as const satisfies Record<string, keyof typeof CORNERS>;

export function tailwindTheme(): string {
  const colors = [...Object.keys(BASE_COLORS), ...Object.keys(ALIAS_COLORS)].map(
    (name) => `--color-${name}: var(--${name});`,
  );

  const fonts = [...Object.keys(TYPEFACES), "heading"].map((name) => `--font-${name}: var(--typeface-${name});`);
  const radii = Object.entries(TAILWIND_RADII).map(([name, corner]) => `--radius-${name}: var(--corner-${corner});`);

  return [
    "/* Generated from packages/design/src/tokens/tokens.ts by `just gen design`; edit the token source instead */",
    "",
    "@theme inline {",
    "  --color-*: initial;",
    "  --color-white: #fff;",
    "  --color-black: #000;",
    ...[...colors, ...fonts, ...radii].map((declaration) => `  ${declaration}`),
    "}",
    "",
  ].join("\n");
}
