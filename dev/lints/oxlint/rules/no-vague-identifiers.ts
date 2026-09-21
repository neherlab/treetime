import { defineRule } from "@oxlint/plugins";
import type { ESTree } from "@oxlint/plugins";

import { declaredIdentifier, VAGUE_NAMES } from "./ast.ts";

export const noVagueIdentifiersRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description:
        "Disallow grab-bag identifiers such as `utils`, `helpers`, `misc`, and `tmp` that name no single responsibility.",
    },
    messages: {
      vague: "`{{name}}` is a grab-bag name. Name the module or binding for its single responsibility.",
    },
  },
  createOnce(context) {
    function check(node: ESTree.Node): void {
      const id = declaredIdentifier(node);

      if (id != null && VAGUE_NAMES.has(id.name)) {
        context.report({ node: id, messageId: "vague", data: { name: id.name } });
      }
    }

    function checkImport(local: ESTree.BindingIdentifier): void {
      if (VAGUE_NAMES.has(local.name)) {
        context.report({ node: local, messageId: "vague", data: { name: local.name } });
      }
    }

    return {
      VariableDeclarator: check,
      FunctionDeclaration: check,
      ClassDeclaration: check,
      TSInterfaceDeclaration: check,
      TSTypeAliasDeclaration: check,
      TSModuleDeclaration: check,
      ImportDefaultSpecifier(node) {
        checkImport(node.local);
      },
      ImportNamespaceSpecifier(node) {
        checkImport(node.local);
      },
    };
  },
});
