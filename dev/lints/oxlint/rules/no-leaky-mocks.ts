import { defineRule } from "@oxlint/plugins";

import { memberChain } from "./ast.ts";

export const noLeakyMocksRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow bare `vi.spyOn` and `vi.fn` that leak into the next test.",
    },
    messages: {
      leakyMock: "Bare vi.spyOn/vi.fn leaks into the next test. Register cleanup or pass an explicit stub.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        const chain = memberChain(node);

        if (chain[0] === "vi" && (chain.at(-1) === "spyOn" || chain.at(-1) === "fn")) {
          context.report({ node, messageId: "leakyMock" });
        }
      },
    };
  },
});
