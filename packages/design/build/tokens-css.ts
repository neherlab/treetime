import { cssColor } from "../src/tokens/color";
import { ALIAS_COLORS, BASE_COLORS, CORNERS, MOTION, TYPEFACES } from "../src/tokens/tokens";

export function tokensCss(): string {
  const declarations = [
    ...Object.entries(BASE_COLORS).map(
      ([name, { light, dark }]) => `--${name}: light-dark(${cssColor(light)}, ${cssColor(dark)});`,
    ),
    ...Object.entries(ALIAS_COLORS).map(([name, target]) => `--${name}: var(--${target});`),
    ...Object.entries(TYPEFACES).map(([name, stack]) => `--typeface-${name}: ${stack.join(", ")};`),
    "--typeface-heading: var(--typeface-sans);",
    ...Object.entries(CORNERS).map(([name, radius]) => `--corner-${name}: ${radius};`),
    `--motion-duration: ${MOTION.duration};`,
    `--motion-easing: ${MOTION.easing};`,
  ];

  return [
    "/* Generated from packages/design/src/tokens/tokens.ts by `just gen design`; edit the token source instead */",
    "",
    ":root {",
    "  color-scheme: light dark;",
    "}",
    "",
    ':root.light,\n:root[data-theme="light"] {',
    "  color-scheme: light;",
    "}",
    "",
    ':root.dark,\n:root[data-theme="dark"] {',
    "  color-scheme: dark;",
    "}",
    "",
    ":root {",
    ...declarations.map((declaration) => `  ${declaration}`),
    "}",
    "",
    "@media (prefers-reduced-motion: reduce) {",
    "  :root {",
    "    --motion-duration: 0ms;",
    "  }",
    "}",
    "",
  ].join("\n");
}
