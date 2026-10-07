import { defineRule } from "@oxlint/plugins";
import * as z from "zod";

const CN_SOURCES = new Set(["clsx", "tailwind-merge", "cn"]);

const zOptions = z.object({ module: z.string() });

export const useThemedCnRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Require importing `cn` from the themed ui module, not `clsx`, `tailwind-merge`, or a bare `cn`.",
    },
    messages: {
      themedCn: "Import `cn` from the themed ui/cn module, not `{{source}}` directly.",
    },
    schema: [
      {
        type: "object",
        properties: { module: { type: "string" } },
        required: ["module"],
        additionalProperties: false,
      },
    ],
  },
  createOnce(context) {
    return {
      ImportDeclaration(node) {
        const options = zOptions.safeParse(context.options[0]);
        const isThemedModule = options.success && context.filename.replaceAll("\\", "/").endsWith(options.data.module);

        if (isThemedModule) {
          return;
        }

        if (typeof node.source.value === "string" && CN_SOURCES.has(node.source.value)) {
          context.report({ node, messageId: "themedCn", data: { source: node.source.value } });
        }
      },
    };
  },
});
