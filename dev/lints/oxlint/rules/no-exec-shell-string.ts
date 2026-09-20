import { defineRule } from "@oxlint/plugins"

import { calleeName } from "./ast.ts"

export const noExecShellStringRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow `exec` and `execSync` with a shell string argument.",
    },
    messages: {
      execString:
        "exec with a shell string is banned. Use execFile or spawn with an argument array.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        const name = calleeName(node)
        if ((name === "exec" || name === "execSync") && node.arguments.length > 0) {
          const first = node.arguments[0]
          if (first.type === "Literal" || first.type === "TemplateLiteral") {
            context.report({ node, messageId: "execString" })
          }
        }
      },
    }
  },
})
