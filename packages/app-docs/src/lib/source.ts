import { type CollectionEntry, getCollection } from "astro:content";

import { CONTENT_ROOT, createDocsSource } from "./docs-source";

const entries = await getCollection("docs");

const metas = await getCollection("meta");

export const entriesById = new Map(entries.map((entry) => [entry.id, entry]));

export const source = createDocsSource(
  {
    pages: entries.map((entry) => ({
      path: contentPath(entry),
      data: { ...entry.data, entryId: entry.id, file: `${CONTENT_ROOT}${contentPath(entry)}` },
    })),
    metas: metas.map((meta) => ({ path: contentPath(meta), data: meta.data })),
  },
  import.meta.env.BASE_URL,
);

function contentPath(entry: CollectionEntry<"docs"> | CollectionEntry<"meta">): string {
  const path = entry.filePath ?? "";

  if (!path.startsWith(CONTENT_ROOT)) {
    throw new Error(`the docs entry ${entry.id} lies outside ${CONTENT_ROOT}: ${path}`);
  }

  return path.slice(CONTENT_ROOT.length);
}
