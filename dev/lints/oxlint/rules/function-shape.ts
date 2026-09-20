import type { ESTree, SourceCode } from "@oxlint/plugins"

export type FunctionNode = ESTree.ArrowFunctionExpression | ESTree.Function

export function isFunctionNode(node: ESTree.Node): node is FunctionNode {
  return (
    node.type === "ArrowFunctionExpression" ||
    node.type === "FunctionDeclaration" ||
    node.type === "FunctionExpression"
  )
}

export function enclosingFunction(node: ESTree.Node): FunctionNode | undefined {
  let current: ESTree.Node | null = node.parent
  while (current !== null) {
    if (isFunctionNode(current)) {
      return current
    }
    current = current.parent
  }
  return undefined
}

export function functionName(node: FunctionNode, sourceCode: SourceCode): string | undefined {
  if (node.id !== null) {
    return node.id.name
  }
  const { parent } = node
  if (parent.type === "VariableDeclarator" && parent.id.type === "Identifier") {
    return parent.id.name
  }
  if (
    parent.type === "MethodDefinition" ||
    parent.type === "PropertyDefinition" ||
    parent.type === "Property"
  ) {
    return keyName(parent.key, sourceCode)
  }
  return undefined
}

export function memberRoot(expression: ESTree.Expression): ESTree.Expression {
  let current = expression
  while (current.type === "MemberExpression" || current.type === "ChainExpression") {
    current = current.type === "ChainExpression" ? current.expression : current.object
  }
  return current
}

function keyName(key: ESTree.PropertyKey, sourceCode: SourceCode): string {
  if (key.type === "Identifier" || key.type === "PrivateIdentifier") {
    return key.name
  }
  if (key.type === "Literal") {
    return String(key.value)
  }
  return sourceCode.getText(key)
}
