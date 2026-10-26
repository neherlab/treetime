import { readFile } from "node:fs/promises";
import { join } from "node:path";
import { fileURLToPath } from "node:url";

import { formatEntry, translationStatus } from "./translation-status";

const DOCS = fileURLToPath(new URL("../../../../docs/", import.meta.url));

const paths = await Array.fromAsync(new Bun.Glob("**/*.{md,mdx,json}").scan({ cwd: DOCS }));

const files = new Map(
  await Promise.all(paths.map(async (path) => [path, await readFile(join(DOCS, path), "utf8")] as const)),
);

for (const entry of translationStatus(files)) {
  process.stdout.write(`${formatEntry(entry)}\n`);
}
