import type { JsonValue, SettingSpec } from "@neherlab/app-contracts";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@neherlab/design/ui/collapsible";
import type { ReactNode } from "react";
import Markdown, { type Components } from "react-markdown";
import rehypeSanitize from "rehype-sanitize";
import remarkGfm from "remark-gfm";

import { m, useLocale, type Locale } from "../i18n/locale";
import { documentedDefault } from "./config";
import { isJsonObject } from "./json";
import { formatList } from "./lists";

const REMARK_PLUGINS = [remarkGfm, remarkHtmlAsText];

const REHYPE_PLUGINS = [rehypeSanitize];

const CODE_CLASS = "bg-muted text-foreground rounded-sm px-1 py-px font-mono text-[0.92em]";

const INLINE_COMPONENTS: Components = {
  p: ({ children }) => <span>{children}</span>,
  code: ({ children }) => <code className={CODE_CLASS}>{children}</code>,
};

const BLOCK_COMPONENTS: Components = {
  p: ({ children }) => <p className="mt-1.5">{children}</p>,
  code: ({ children }) => <code className={CODE_CLASS}>{children}</code>,
};

export function SettingSummary({ spec, label }: { spec: SettingSpec; label: ReactNode }) {
  const locale = useLocale();

  return (
    <div className="flex flex-wrap items-baseline gap-2">
      {label} <code className="text-muted-foreground font-mono text-xs">{spec.flag}</code>{" "}
      {spec.role === "setting" && (
        <span className="text-muted-foreground text-xs">
          {m.settings_default({ value: defaultText(documentedDefault(spec), locale) }, { locale })}
        </span>
      )}
    </div>
  );
}

export function SettingHelp({ spec }: { spec: SettingSpec }) {
  const locale = useLocale();

  return (
    <div className="text-muted-foreground max-w-[72ch] text-[0.8125rem]">
      <Markdown remarkPlugins={REMARK_PLUGINS} rehypePlugins={REHYPE_PLUGINS} components={INLINE_COMPONENTS}>
        {spec.help}
      </Markdown>
      {spec.more !== "" && (
        <Collapsible className="inline">
          <CollapsibleTrigger className="text-primary ml-1 underline-offset-4 hover:underline">
            {m.settings_more({}, { locale })}
          </CollapsibleTrigger>
          <CollapsibleContent hiddenUntilFound>
            <Markdown remarkPlugins={REMARK_PLUGINS} rehypePlugins={REHYPE_PLUGINS} components={BLOCK_COMPONENTS}>
              {spec.more}
            </Markdown>
          </CollapsibleContent>
        </Collapsible>
      )}
    </div>
  );
}

export function defaultText(value: JsonValue | undefined, locale: Locale): string {
  if (value === undefined) {
    return m.settings_not_set({}, { locale });
  }

  if (Array.isArray(value)) {
    return value.length === 0 ? m.settings_none({}, { locale }) : formatList(value);
  }

  if (value === true) {
    return m.settings_on({}, { locale });
  }

  if (value === false) {
    return m.settings_off({}, { locale });
  }

  return isJsonObject(value) ? JSON.stringify(value) : String(value);
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
