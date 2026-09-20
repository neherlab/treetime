import { defineRule } from "@oxlint/plugins"

import { rootCalleeName, TEST_CALLERS } from "./ast.ts"

export const noFakeSuccessRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow a test whose callback body is empty and therefore always passes.",
    },
    messages: {
      empty: "A test with an empty body always passes. Assert observable behavior or delete it.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        if (!TEST_CALLERS.has(rootCalleeName(node) ?? "")) {
          return
        }
        const body = node.arguments.at(-1)
        if (
          body != null &&
          (body.type === "ArrowFunctionExpression" || body.type === "FunctionExpression") &&
          body.body.type === "BlockStatement" &&
          body.body.body.length === 0
        ) {
          context.report({ node, messageId: "empty" })
        }
      },
    }
  },
})
