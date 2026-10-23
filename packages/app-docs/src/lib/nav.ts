import type * as PageTree from "fumadocs-core/page-tree";
import type { ReactNode } from "react";
import { match } from "ts-pattern";

export function navNodes(nodes: readonly PageTree.Node[]): NavNode[] {
  return nodes.map((node) =>
    match(node)
      .returnType<NavNode>()
      .with({ type: "page" }, (page) => ({ kind: "page", name: nameText(page.name), url: page.url }))
      .with({ type: "folder" }, (folder) => ({
        kind: "folder",
        name: nameText(folder.name),
        url: folder.index?.url,
        children: navNodes(folder.children),
      }))
      .with({ type: "separator" }, (separator) => ({ kind: "separator", name: nameText(separator.name ?? "") }))
      .exhaustive(),
  );
}

export function containsUrl(node: NavNode, url: string): boolean {
  return match(node)
    .with({ kind: "page" }, (page) => page.url === url)
    .with(
      { kind: "folder" },
      (folder) => folder.url === url || folder.children.some((child) => containsUrl(child, url)),
    )
    .with({ kind: "separator" }, () => false)
    .exhaustive();
}

export type NavNode =
  | { kind: "page"; name: string; url: string }
  | { kind: "folder"; name: string; url: string | undefined; children: NavNode[] }
  | { kind: "separator"; name: string };

function nameText(name: ReactNode): string {
  // oxlint-disable-next-line anti-slop/no-runtime-typeof -- the page tree types names as React nodes, and the loader fills them from title strings
  if (typeof name !== "string") {
    throw new TypeError("a docs navigation entry has a name that is no text");
  }

  return name;
}
