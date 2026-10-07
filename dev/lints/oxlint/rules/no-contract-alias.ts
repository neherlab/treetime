import { defineRule } from "@oxlint/plugins";
import type { ESTree } from "@oxlint/plugins";

import { contractImports, contractsPackage, isContractsSource } from "./contracts.ts";

export const noContractAliasRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description:
        "Disallow aliases and patched copies of generated contract types; import the generated name directly.",
    },
    messages: {
      alias: "`{{alias}}` only renames the generated type `{{name}}`. Import `{{name}}` directly.",
      patched:
        "`Omit<{{name}}, ...> & {...}` re-declares part of the generated type `{{name}}`. Change the Rust type, then regenerate.",
      renamedExport: "Re-exporting the generated `{{name}}` as `{{alias}}` gives one type two names. Export it as is.",
    },
    schema: [
      {
        type: "object",
        properties: { package: { type: "string" } },
        required: ["package"],
        additionalProperties: false,
      },
    ],
  },
  createOnce(context) {
    let contracts = new Set<string>();
    let contractsSource: string | undefined;

    return {
      Program(node) {
        contractsSource = contractsPackage(context.options);
        contracts = contractsSource === undefined ? new Set() : contractImports(node, contractsSource);
      },
      TSTypeAliasDeclaration(node) {
        const name = referencedName(node.typeAnnotation);

        if (name !== undefined && contracts.has(name)) {
          context.report({ node, messageId: "alias", data: { alias: node.id.name, name } });
        }
      },
      TSIntersectionType(node) {
        const omitted = node.types.flatMap((member) => {
          const name = omittedName(member);

          return name !== undefined && contracts.has(name) ? [name] : [];
        });

        const [name] = omitted;

        if (name !== undefined && node.types.some((member) => member.type === "TSTypeLiteral")) {
          context.report({ node, messageId: "patched", data: { name } });
        }
      },
      ExportNamedDeclaration(node) {
        const fromContracts =
          node.source !== null &&
          contractsSource !== undefined &&
          isContractsSource(node.source.value, contractsSource);

        for (const specifier of node.specifiers) {
          const local = exportName(specifier.local);
          const exported = exportName(specifier.exported);

          if (local !== exported && (fromContracts || contracts.has(local))) {
            context.report({ node: specifier, messageId: "renamedExport", data: { name: local, alias: exported } });
          }
        }
      },
    };
  },
});

function referencedName(type: ESTree.TSType): string | undefined {
  return type.type === "TSTypeReference" && type.typeName.type === "Identifier" && type.typeArguments == null
    ? type.typeName.name
    : undefined;
}

function omittedName(type: ESTree.TSType): string | undefined {
  if (type.type !== "TSTypeReference" || type.typeName.type !== "Identifier" || type.typeName.name !== "Omit") {
    return undefined;
  }

  const [target] = type.typeArguments?.params ?? [];

  return target?.type === "TSTypeReference" && target.typeName.type === "Identifier" ? target.typeName.name : undefined;
}

function exportName(node: ESTree.ModuleExportName): string {
  return node.type === "Literal" ? node.value : node.name;
}
