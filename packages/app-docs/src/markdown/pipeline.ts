import { unified } from "@astrojs/markdown-remark";
import { remarkDirectiveAdmonition } from "fumadocs-core/mdx-plugins/remark-directive-admonition";
import { remarkSteps } from "fumadocs-core/mdx-plugins/remark-steps";
import rehypeKatex from "rehype-katex";
import remarkDirective from "remark-directive";
import remarkMath from "remark-math";

import { TEX_MACROS } from "../math/macros";
import { rehypeBaseUrls } from "./rehype-base-urls";
import { remarkProvidedComponents } from "./remark-provided-components";
import { restoreDirectives } from "./restore-directives";

const CALLOUT_TYPES = { note: "info", tip: "success", warning: "warning", danger: "error" };

export function docsMarkdownProcessor(base: string, contentRoots: readonly string[]) {
  return unified({
    remarkPlugins: [
      remarkDirective,
      [remarkDirectiveAdmonition, { types: CALLOUT_TYPES }],
      restoreDirectives,
      remarkSteps,
      remarkMath,
      remarkProvidedComponents,
    ],
    rehypePlugins: [
      [rehypeKatex, { macros: TEX_MACROS }],
      [rehypeBaseUrls, { base, contentRoots }],
    ],
  });
}
