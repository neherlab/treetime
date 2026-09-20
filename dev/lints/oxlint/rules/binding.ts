import type { ESTree, Scope, SourceCode, Variable } from "@oxlint/plugins";

export function bindingVariable(identifier: ESTree.IdentifierReference, sourceCode: SourceCode): Variable | undefined {
  let scope: Scope | null = sourceCode.getScope(identifier);
  while (scope !== null) {
    const reference = scope.references.find((entry) => entry.identifier === identifier);
    if (reference !== undefined) {
      return reference.resolved ?? undefined;
    }
    scope = scope.upper;
  }
  return undefined;
}

export function bindingImport(
  identifier: ESTree.IdentifierReference,
  sourceCode: SourceCode,
): ImportBinding | undefined {
  const definition = bindingVariable(identifier, sourceCode)?.defs[0];
  if (definition?.type !== "ImportBinding") {
    return undefined;
  }
  const { node } = definition;
  if (
    node.type !== "ImportSpecifier" &&
    node.type !== "ImportNamespaceSpecifier" &&
    node.type !== "ImportDefaultSpecifier"
  ) {
    return undefined;
  }
  const declaration = node.parent;
  if (
    declaration.type !== "ImportDeclaration" ||
    declaration.importKind === "type" ||
    (node.type === "ImportSpecifier" && node.importKind === "type")
  ) {
    return undefined;
  }
  return { node, source: declaration.source.value };
}

export interface ImportBinding {
  node: ESTree.ImportSpecifier | ESTree.ImportNamespaceSpecifier | ESTree.ImportDefaultSpecifier;
  source: string;
}

export function bindingImportName(binding: ImportBinding): string | undefined {
  if (binding.node.type === "ImportDefaultSpecifier") {
    return "default";
  }
  if (binding.node.type === "ImportSpecifier") {
    const { imported } = binding.node;
    return imported.type === "Identifier" ? imported.name : imported.value;
  }
  return undefined;
}

export function bindingInitializer(
  identifier: ESTree.IdentifierReference,
  sourceCode: SourceCode,
): ESTree.Expression | undefined {
  const variable = bindingVariable(identifier, sourceCode);
  const definition = variable?.defs[0];
  if (
    definition?.type !== "Variable" ||
    definition.node.type !== "VariableDeclarator" ||
    definition.node.id.type !== "Identifier" ||
    definition.node.end >= identifier.start ||
    definition.parent?.type !== "VariableDeclaration" ||
    definition.parent.kind === "var" ||
    variable?.references.some((reference) => reference.isWrite() && !reference.init) === true
  ) {
    return undefined;
  }
  return definition.node.init ?? undefined;
}
