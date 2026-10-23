import { locales, type Locale } from "@neherlab/app-ui/i18n";
import type * as PageTree from "fumadocs-core/page-tree";
import { findNeighbour } from "fumadocs-core/page-tree";

import type { DocsPageData } from "./docs-source";
import { isFallbackPage } from "./docs-source";
import { navNodes, type NavNode } from "./nav";
import { docsPath } from "./urls";

const FOOTNOTE_LABEL = "footnote-label";

const EDIT_ROOT = "https://github.com/neherlab/treetime/edit/rust/packages/app-docs/";

export function pageView({ page, locale, tree, headings, base }: PageViewInput): PageView {
  const neighbours = findNeighbour(tree, page.url);

  const tocHeadings = headings.filter(
    (heading) => (heading.depth === 2 || heading.depth === 3) && heading.slug !== FOOTNOTE_LABEL,
  );

  return {
    nodes: navNodes(tree.children),
    currentUrl: page.url,
    version: base.replaceAll(/^\/+|\/+$/gu, ""),
    slug: page.slugs.join("/"),
    localeUrls: new Map(locales.map((target) => [target, docsPath(base, page.slugs, target)])),
    tocHeadings,
    showToc: page.data.tableOfContents && tocHeadings.length > 0,
    fallback: isFallbackPage(page, locale),
    previous: footerLink(neighbours.previous),
    next: footerLink(neighbours.next),
    editUrl: page.data.editUrl === false ? undefined : `${EDIT_ROOT}${page.data.file}`,
  };
}

export interface PageViewInput {
  page: { url: string; path: string; slugs: string[]; data: DocsPageData };
  locale: Locale;
  tree: PageTree.Root;
  headings: ReadonlyArray<{ depth: number; slug: string; text: string }>;
  base: string;
}

export interface PageView {
  nodes: NavNode[];
  currentUrl: string;
  version: string;
  slug: string;
  localeUrls: ReadonlyMap<Locale, string>;
  tocHeadings: Array<{ depth: number; slug: string; text: string }>;
  showToc: boolean;
  fallback: boolean;
  previous: { name: string; url: string } | undefined;
  next: { name: string; url: string } | undefined;
  editUrl: string | undefined;
}

function footerLink(item: PageTree.Item | undefined): { name: string; url: string } | undefined {
  return item === undefined ? undefined : { name: navNodes([item])[0]?.name ?? item.url, url: item.url };
}
