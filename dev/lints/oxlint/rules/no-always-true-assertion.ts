import { defineRule } from "@oxlint/plugins"

import { calleeName } from "./ast.ts"

export const noAlwaysTrueAssertionRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow `expect` on a literal, which asserts an always-true value.",
    },
    messages: {
      literal: "Asserting on a literal is always true. Assert on a computed value.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        if (calleeName(node) !== "expect" || node.arguments.length === 0) {
          return
        }
        if (node.arguments[0].type === "Literal") {
          context.report({ node, messageId: "literal" })
        }
      },
    }
  },
})
