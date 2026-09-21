import { defineRule } from "@oxlint/plugins";

import { TYPOGRAPHIC } from "./ast.ts";

export const noTypographicCharactersRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Disallow curly quotes, en and em dashes, and emoji in string and template text.",
    },
    messages: {
      typographic: "Curly quotes, dashes, and emoji are banned. Use straight ASCII punctuation.",
    },
  },
  createOnce(context) {
    return {
      Literal(node) {
        if (typeof node.value === "string" && TYPOGRAPHIC.test(node.value)) {
          context.report({ node, messageId: "typographic" });
        }
      },
      TemplateElement(node) {
        if (typeof node.value.raw === "string" && TYPOGRAPHIC.test(node.value.raw)) {
          context.report({ node, messageId: "typographic" });
        }
      },
    };
  },
});
