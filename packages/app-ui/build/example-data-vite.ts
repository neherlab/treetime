import { existsSync, readFileSync } from "node:fs";
import { join, posix } from "node:path";
import { fileURLToPath } from "node:url";

import { type DataIndex, zDataIndex } from "@neherlab/app-contracts";
import { convertPathToPattern } from "tinyglobby";
import type { Plugin } from "vite";
import { viteStaticCopy } from "vite-plugin-static-copy";

const EXAMPLE_DATA_FOLDER = "data";

const DATA_INDEX_FILE = "index.json";

const DATA_DIR = fileURLToPath(new URL("../../../data", import.meta.url));

export function exampleData(dir: string = DATA_DIR): Plugin[] {
  const index = parseDataIndex(readFileSync(join(dir, DATA_INDEX_FILE), "utf8"), (path) => existsSync(join(dir, path)));

  return viteStaticCopy({ targets: exampleDataTargets(dir, index) });
}

export function parseDataIndex(text: string, exists: (path: string) => boolean): DataIndex {
  const parsed = zDataIndex.safeParse(JSON.parse(text));

  if (!parsed.success) {
    throw new Error(`The example data index does not match the DataIndex schema: ${parsed.error.message}`);
  }

  const paths = parsed.data.files.map((file) => file.path);
  const outside = paths.filter((path) => !isDataPath(path));

  if (outside.length > 0) {
    throw new Error(`The example data index lists paths outside the data folder: ${outside.join(", ")}`);
  }

  const missing = paths.filter((path) => !exists(path));

  if (missing.length > 0) {
    throw new Error(
      `The example data index lists missing files: ${missing.join(", ")}. Regenerate it with: just data-reindex`,
    );
  }

  return parsed.data;
}

export function exampleDataTargets(dir: string, index: DataIndex) {
  return [DATA_INDEX_FILE, ...index.files.map((file) => file.path)].map((path) => ({
    src: convertPathToPattern(join(dir, path)),
    dest: posix.join(EXAMPLE_DATA_FOLDER, posix.dirname(path)),
    rename: { stripBase: true } as const,
  }));
}

function isDataPath(path: string): boolean {
  return path !== "" && posix.normalize(path) === path && !posix.isAbsolute(path) && !path.startsWith("../");
}
