import { defineRule } from "@oxlint/plugins"

import { memberChain } from "./ast.ts"

export const noModuleMocksRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow `vi.mock`, `vi.doMock`, and `vi.unmock` module mocking.",
    },
    messages: {
      moduleMock: "`{{chain}}` module mocking is banned. Inject collaborators instead.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        const chain = memberChain(node)
        const method = chain.at(-1)
        if (chain[0] === "vi" && (method === "mock" || method === "doMock" || method === "unmock")) {
          context.report({ node, messageId: "moduleMock", data: { chain: chain.join(".") } })
        }
      },
    }
  },
})
