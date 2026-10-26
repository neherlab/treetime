import type { Element, Root } from "hast";
import type { MdxJsxFlowElementHast, MdxJsxTextElementHast } from "mdast-util-mdx-jsx";
import { visit } from "unist-util-visit";
import type { VFile } from "vfile";

import { assetPath, docsPath, docsLocale } from "../lib/urls";

const URL_PROPERTIES = ["href", "src"] as const;

const PAGE_FILE = /\.(?:md|mdx)$/u;

const PUBLIC_FILE = /\.[a-z0-9]+$/iu;

export function rehypeBaseUrls({ base, contentRoots }: { base: string; contentRoots: readonly string[] }) {
  return (tree: Root, file: VFile): void => {
    const locale = pageLocale(file.path, contentRoots);

    visit(tree, (node) => {
      if (node.type === "element") {
        rewriteProperties(node, base, locale);
      } else if (node.type === "mdxJsxFlowElement" || node.type === "mdxJsxTextElement") {
        rewriteAttributes(node, base, locale);
      }
    });
  };
}

export function withBase(url: string, base: string, locale: string | undefined): string {
  const { pathname, rest } = splitUrl(url);

  if (PUBLIC_FILE.test(pathname) && !PAGE_FILE.test(pathname)) {
    return `${assetPath(base, pathname)}${rest}`;
  }

  return `${docsPath(
    base,
    pathname.split("/").filter((segment) => segment !== ""),
    docsLocale(locale),
  )}${rest}`;
}

function rewriteProperties(node: Element, base: string, locale: string | undefined): void {
  for (const property of URL_PROPERTIES) {
    const value = node.properties[property];

    // oxlint-disable-next-line anti-slop/no-runtime-typeof -- hast property values are strings, numbers, booleans, or lists, and only strings hold URLs
    if (typeof value === "string" && isRootRelative(value)) {
      node.properties[property] = withBase(value, base, locale);
    }
  }
}

function rewriteAttributes(
  node: MdxJsxFlowElementHast | MdxJsxTextElementHast,
  base: string,
  locale: string | undefined,
): void {
  for (const attribute of node.attributes) {
    if (
      attribute.type === "mdxJsxAttribute" &&
      URL_PROPERTIES.some((property) => property === attribute.name) &&
      // oxlint-disable-next-line anti-slop/no-runtime-typeof -- MDX attribute values are strings or expressions, and only strings hold URLs
      typeof attribute.value === "string" &&
      isRootRelative(attribute.value)
    ) {
      attribute.value = withBase(attribute.value, base, locale);
    }
  }
}

function isRootRelative(url: string): boolean {
  return url.startsWith("/") && !url.startsWith("//");
}

function splitUrl(url: string): { pathname: string; rest: string } {
  const end = url.search(/[?#]/u);

  return end === -1 ? { pathname: url, rest: "" } : { pathname: url.slice(0, end), rest: url.slice(end) };
}

function pageLocale(path: string, contentRoots: readonly string[]): string | undefined {
  const root = contentRoots.find((candidate) => path.startsWith(candidate));

  return root === undefined ? undefined : path.slice(root.length).split("/")[0];
}
