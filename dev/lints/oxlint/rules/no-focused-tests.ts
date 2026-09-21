import { defineRule } from "@oxlint/plugins";

import { memberChain, rootCalleeName, TEST_CALLERS } from "./ast.ts";

export const noFocusedTestsRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow focused tests: `fit`, `fdescribe`, and `.only`.",
    },
    messages: {
      focused: "Focused tests are banned. They hide the rest of the suite.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        const root = rootCalleeName(node);

        if (root === "fit" || root === "fdescribe") {
          context.report({ node, messageId: "focused" });

          return;
        }

        if (TEST_CALLERS.has(root ?? "") && memberChain(node).includes("only")) {
          context.report({ node, messageId: "focused" });
        }
      },
    };
  },
});
