import { readFile } from "node:fs/promises";
import { join } from "node:path";

import { close, createIndex } from "pagefind";

import { pageUrl } from "./built-pages";

export async function writeSearchIndex(root: string, base: string, files: readonly string[]): Promise<void> {
  try {
    const created = await createIndex({});
    failOnErrors("create the search index", created.errors);
    const { index } = created;

    if (index === undefined) {
      throw new Error("Pagefind created no search index");
    }

    await Promise.all(
      files.map(async (file) => {
        const content = await readFile(join(root, file), "utf8");
        const added = await index.addHTMLFile({ url: pageUrl(base, file), content });
        failOnErrors(`index ${file}`, added.errors);
      }),
    );

    const written = await index.writeFiles({ outputPath: join(root, "pagefind") });
    failOnErrors("write the search index", written.errors);
  } finally {
    await close();
  }
}

function failOnErrors(step: string, errors: readonly string[]): void {
  if (errors.length > 0) {
    throw new Error(`Pagefind could not ${step}: ${errors.join("; ")}`);
  }
}
