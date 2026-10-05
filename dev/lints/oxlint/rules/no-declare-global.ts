import { defineRule } from "@oxlint/plugins";

export const noDeclareGlobalRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow `declare global` blocks, which widen the global types of every module in the program.",
    },
    messages: {
      declareGlobal:
        "`declare global` changes the global types of every module. Pass the value through a typed module or context instead.",
    },
  },
  createOnce(context) {
    return {
      TSModuleDeclaration(node) {
        if (node.kind === "global") {
          context.report({ node, messageId: "declareGlobal" });
        }
      },
    };
  },
});
