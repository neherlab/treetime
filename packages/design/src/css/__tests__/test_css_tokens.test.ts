import { ALIAS_COLORS, BASE_COLORS, CORNERS, TYPEFACES } from "../../tokens/tokens";

import docusaurusCss from "../docusaurus.css?raw";
import revealCss from "../reveal.css?raw";
import tailwindCss from "../tailwind.css?raw";

const TOKEN_NAMES = new Set([
  ...Object.keys(BASE_COLORS),
  ...Object.keys(ALIAS_COLORS),
  ...[...Object.keys(TYPEFACES), "heading"].map((name) => `typeface-${name}`),
  ...Object.keys(CORNERS).map((name) => `corner-${name}`),
  "motion-duration",
  "motion-easing",
]);

const FOREIGN_PREFIXES = ["ifm-", "r-", "tw-", "docusaurus-"];

describe("design stylesheets", () => {
  test.each([
    ["docusaurus.css", docusaurusCss],
    ["reveal.css", revealCss],
    ["tailwind.css", tailwindCss],
  ])("%s reads only variables that the tokens define", (_name, css) => {
    const unknown = [...css.matchAll(/var\(--(?<name>[a-z0-9-]+)/gu)].flatMap((match) => {
      const name = match.groups?.["name"] ?? "";

      return TOKEN_NAMES.has(name) || FOREIGN_PREFIXES.some((prefix) => name.startsWith(prefix)) ? [] : [name];
    });

    expect(unknown).toStrictEqual([]);
  });
});
