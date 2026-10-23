import type { Parent, Root } from "mdast";
import type { ContainerDirective, LeafDirective, TextDirective } from "mdast-util-directive";
import { directiveToMarkdown } from "mdast-util-directive";
import { toMarkdown } from "mdast-util-to-markdown";
import { SKIP, visit } from "unist-util-visit";
import type { VFile } from "vfile";

export function restoreDirectives() {
  return (tree: Root, file: VFile): void => {
    visit(tree, (node, index, parent) => {
      if (index === undefined || parent === undefined || !isDirective(node) || node.data !== undefined) {
        return undefined;
      }

      if (node.type === "containerDirective") {
        file.fail(`The block \`:::${node.name}\` is no callout type. Use note, tip, warning, or danger.`, {
          place: node.position,
        });
      }

      parent.children[index] = restoredNode(node);

      return SKIP;
    });
  };
}

function isDirective(node: { type: string }): node is TextDirective | LeafDirective | ContainerDirective {
  return node.type === "textDirective" || node.type === "leafDirective" || node.type === "containerDirective";
}

function restoredNode(node: TextDirective | LeafDirective | ContainerDirective): Parent["children"][number] {
  const text = {
    type: "text" as const,
    value: toMarkdown(node, { extensions: [directiveToMarkdown()] }).replace(/\n$/u, ""),
  };

  return node.type === "textDirective" ? text : { type: "paragraph", children: [text] };
}
