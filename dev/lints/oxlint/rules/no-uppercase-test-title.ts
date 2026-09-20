import { defineRule } from "@oxlint/plugins"

import { rootCalleeName, TEST_CALLERS } from "./ast.ts"

export const noUppercaseTestTitleRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Require test titles to start lowercase and read as a sentence fragment.",
    },
    messages: {
      uppercase: "Test titles must start lowercase and read as a sentence fragment.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        if (!TEST_CALLERS.has(rootCalleeName(node) ?? "") || node.arguments.length === 0) {
          return
        }
        const title = node.arguments[0]
        if (title.type === "Literal" && typeof title.value === "string" && /^[A-Z]/.test(title.value)) {
          context.report({ node: title, messageId: "uppercase" })
        }
      },
    }
  },
})
