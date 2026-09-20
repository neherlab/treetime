import { defineRule } from "@oxlint/plugins"

export const useClassNameHelperRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Disallow importing `clsx` or `tailwind-merge` directly; use the themed cn helper.",
    },
    messages: {
      useCn: "Import the themed cn helper instead of `{{source}}` directly.",
    },
  },
  createOnce(context) {
    return {
      ImportDeclaration(node) {
        const source = node.source.value
        if (source === "clsx" || source === "tailwind-merge") {
          context.report({ node, messageId: "useCn", data: { source } })
        }
      },
    }
  },
})
