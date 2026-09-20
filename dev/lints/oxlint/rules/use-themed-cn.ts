import { defineRule } from "@oxlint/plugins"

const CN_SOURCES = new Set(["clsx", "tailwind-merge", "cn"])

export const useThemedCnRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Require importing `cn` from the themed ui module, not `clsx`, `tailwind-merge`, or a bare `cn`.",
    },
    messages: {
      themedCn: "Import `cn` from the themed ui/cn module, not `{{source}}` directly.",
    },
  },
  createOnce(context) {
    const filename = context.filename
    const isThemedModule = filename.replace(/\\/g, "/").endsWith("app-ui/src/ui/cn.ts")
    return {
      ImportDeclaration(node) {
        if (isThemedModule) {
          return
        }
        if (typeof node.source.value === "string" && CN_SOURCES.has(node.source.value)) {
          context.report({ node, messageId: "themedCn", data: { source: node.source.value } })
        }
      },
    }
  },
})
