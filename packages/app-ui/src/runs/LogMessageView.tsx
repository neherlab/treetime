import { memo, type ReactNode } from "react";
import Markdown, { type Components } from "react-markdown";
import rehypeSanitize from "rehype-sanitize";
import remarkGfm from "remark-gfm";

const REMARK_PLUGINS = [remarkGfm, remarkHtmlAsText];

const REHYPE_PLUGINS = [rehypeSanitize];

const heading = ({ children }: { children?: ReactNode }) => <span className="font-bold">{children}</span>;

const LOG_COMPONENTS: Components = {
  p: ({ children }) => <span>{children}</span>,
  h1: heading,
  h2: heading,
  h3: heading,
  h4: heading,
  h5: heading,
  h6: heading,
  ul: ({ children }) => <ul className="list-disc pl-5 whitespace-normal">{children}</ul>,
  ol: ({ children }) => <ol className="list-decimal pl-5 whitespace-normal">{children}</ol>,
  table: ({ children }) => (
    <table className="my-1 border-separate border-spacing-x-3 text-right whitespace-normal tabular-nums">
      {children}
    </table>
  ),
  th: ({ children }) => <th className="text-muted-foreground font-normal">{children}</th>,
  pre: ({ children }) => <pre className="bg-muted my-1 overflow-x-auto rounded-sm px-2 py-1">{children}</pre>,
  code: ({ children }) => <code className="bg-muted rounded-sm px-1">{children}</code>,
  a: ({ children, href }) => (
    <a
      href={href}
      target="_blank"
      rel="noopener noreferrer"
      className="text-primary underline-offset-4 hover:underline"
    >
      {children}
    </a>
  ),
};

export const LogMessageView = memo(function LogMessageView({ message }: { message: string }) {
  return (
    <Markdown remarkPlugins={REMARK_PLUGINS} rehypePlugins={REHYPE_PLUGINS} components={LOG_COMPONENTS}>
      {message}
    </Markdown>
  );
});

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
