import { defineRule } from "@oxlint/plugins";

import { memberChain, rootCalleeName, TEST_CALLERS } from "./ast.ts";

export const noDisabledTestsRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow disabled tests: `xit`, `xdescribe`, `xtest`, and `.skip`.",
    },
    messages: {
      disabled: "Disabled tests are banned. Fix or delete the test.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        const root = rootCalleeName(node);

        if (root === "xit" || root === "xdescribe" || root === "xtest") {
          context.report({ node, messageId: "disabled" });

          return;
        }

        if (TEST_CALLERS.has(root ?? "") && memberChain(node).includes("skip")) {
          context.report({ node, messageId: "disabled" });
        }
      },
    };
  },
});
