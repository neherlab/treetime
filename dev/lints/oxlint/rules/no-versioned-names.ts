import { defineRule } from "@oxlint/plugins";
import type { ESTree } from "@oxlint/plugins";

import { declaredIdentifier, VERSIONED_NAME } from "./ast.ts";

export const noVersionedNamesRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Disallow version suffixes such as `_v2`, `_old`, `_new`, and `_tmp` on declared names.",
    },
    messages: {
      versioned: "`{{name}}` carries a version suffix. Names must describe the final state, not edit history.",
    },
  },
  createOnce(context) {
    function check(node: ESTree.Node): void {
      const id = declaredIdentifier(node);

      if (id != null && VERSIONED_NAME.test(id.name)) {
        context.report({ node: id, messageId: "versioned", data: { name: id.name } });
      }
    }

    return {
      VariableDeclarator: check,
      FunctionDeclaration: check,
      ClassDeclaration: check,
      TSInterfaceDeclaration: check,
      TSTypeAliasDeclaration: check,
    };
  },
});
