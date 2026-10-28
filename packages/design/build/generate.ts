import { mkdir, writeFile } from "node:fs/promises";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { parseArgs } from "node:util";

import { matplotlibStyle } from "./matplotlib-style";
import { rustPalette } from "./rust-palette";
import { tailwindTheme } from "./tailwind-theme";
import { tokensCss } from "./tokens-css";

const ROOT = fileURLToPath(new URL("../../..", import.meta.url));

const OUTPUTS: readonly Output[] = [
  { path: "packages/design/src/css/tokens.css", render: tokensCss },
  { path: "packages/design/src/css/tailwind-theme.css", render: tailwindTheme },
  { path: "packages/design/matplotlib/neherlab.mplstyle", render: matplotlibStyle },
  { path: "packages/treetime-charts/src/palette.rs", render: rustPalette },
];

const { values } = parseArgs({ options: { outdir: { type: "string", default: ROOT } } });

for (const output of OUTPUTS) {
  const target = join(values.outdir, output.path);
  await mkdir(dirname(target), { recursive: true });
  await writeFile(target, output.render());
  console.log(output.path);
}

interface Output {
  path: string;
  render: () => string;
}
