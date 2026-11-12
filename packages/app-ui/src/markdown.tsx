import Markdown, { type Components } from "react-markdown";
import rehypeSanitize from "rehype-sanitize";
import remarkGfm from "remark-gfm";

const REMARK_PLUGINS = [remarkGfm, remarkHtmlAsText];

const REHYPE_PLUGINS = [rehypeSanitize];

export type { Components as MarkdownComponents };

export function MarkdownText({ text, components }: { text: string; components: Components }) {
  return (
    <Markdown remarkPlugins={REMARK_PLUGINS} rehypePlugins={REHYPE_PLUGINS} components={components}>
      {text}
    </Markdown>
  );
}

function remarkHtmlAsText() {
  return htmlAsText;
}

function htmlAsText(node: MarkdownNode): void {
  if (node.type === "html") {
    node.type = "text";
  }

  for (const child of node.children ?? []) {
    htmlAsText(child);
  }
}

interface MarkdownNode {
  type: string;
  children?: MarkdownNode[];
}
