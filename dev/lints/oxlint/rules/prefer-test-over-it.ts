import { defineRule } from "@oxlint/plugins";

export const preferTestOverItRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Prefer `test` over `it` as the test declaration.",
    },
    messages: {
      useTest: "Use `test` instead of `it`.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        if (node.callee.type === "Identifier" && node.callee.name === "it") {
          context.report({ node: node.callee, messageId: "useTest" });
        }
      },
    };
  },
});
