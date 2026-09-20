import { defineRule } from "@oxlint/plugins"

const FILE_DISABLE = /^\s*oxlint-disable(?!-next-line|-line)\b/
const ENABLE = /^\s*oxlint-enable\b/

export const noUnboundedSuppressionRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Require a matching oxlint-enable after a file-wide oxlint-disable.",
    },
    messages: {
      unbounded:
        "A file-wide oxlint-disable silences the rest of the file. Close its scope with a matching oxlint-enable, or use oxlint-disable-next-line.",
    },
  },
  createOnce(context) {
    return {
      Program() {
        const comments = context.sourceCode.getAllComments()
        const enables = comments.filter((comment) => ENABLE.test(comment.value))
        for (const comment of comments) {
          if (
            FILE_DISABLE.test(comment.value) &&
            !enables.some((enable) => enable.start > comment.end)
          ) {
            context.report({ loc: comment.loc, messageId: "unbounded" })
          }
        }
      },
    }
  },
})
