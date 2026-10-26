import { baseLocale, m } from "@neherlab/app-ui/i18n";
import { COMMANDS } from "@neherlab/app-ui/settings/catalog";
import { type CollectionEntry, getCollection } from "astro:content";

import { CONTENT_ROOTS, type ContentRoot } from "./content-roots";
import { createDocsSource, type DocsPageFile } from "./docs-source";

const pageEntries = [...(await getCollection("docs")), ...(await getCollection("generated"))];

const metas = await getCollection("meta");

export const entriesByKey = new Map(pageEntries.map((entry) => [entryKey(entry), entry]));

export const source = createDocsSource(
  {
    pages: uniquePaths([...pageEntries.map(entryPage), ...settingsPages()]),
    metas: metas.map((meta) => ({ path: contentPath(meta, "docs"), data: meta.data })),
  },
  import.meta.env.BASE_URL,
);

function entryPage(entry: PageEntry): DocsPageFile {
  const path = contentPath(entry, entry.collection);
  const { title, description, template, tableOfContents } = entry.data;

  return {
    path,
    data: {
      title,
      description,
      template,
      tableOfContents,
      body: { kind: "entry", key: entryKey(entry) },
      editPath: entry.collection === "docs" ? path : undefined,
    },
  };
}

function settingsPages(): DocsPageFile[] {
  return COMMANDS.map(({ command, title }) => ({
    path: `${baseLocale}/reference/settings/${command}.mdx`,
    data: {
      title: command,
      description: m.docs_settings_description({ command, title }, { locale: baseLocale }),
      template: "doc",
      tableOfContents: false,
      body: { kind: "settings", command },
      editPath: undefined,
    },
  }));
}

function uniquePaths(pages: DocsPageFile[]): DocsPageFile[] {
  const repeated = pages.flatMap((page, index) =>
    pages.findIndex((other) => other.path === page.path) === index ? [] : [page.path],
  );

  if (repeated.length > 0) {
    throw new Error(`several docs pages have the same path: ${[...new Set(repeated)].join(", ")}`);
  }

  return pages;
}

function entryKey(entry: PageEntry): string {
  return `${entry.collection}:${entry.id}`;
}

function contentPath(entry: PageEntry | CollectionEntry<"meta">, root: ContentRoot): string {
  const path = entry.filePath ?? "";
  const prefix = CONTENT_ROOTS[root];

  if (!path.startsWith(prefix)) {
    throw new Error(`the docs entry ${entry.id} lies outside ${prefix}: ${path}`);
  }

  return path.slice(prefix.length);
}

type PageEntry = CollectionEntry<"docs"> | CollectionEntry<"generated">;
