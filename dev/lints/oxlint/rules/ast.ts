import type { ESTree, Ranged } from "@oxlint/plugins";

export const VERSIONED_NAME = /(?:_v\d+|_old|_new|_fixed|_deprecated|_copy\d*|_backup|_bak|_temp|_tmp)$/iu;

export const VAGUE_NAMES = new Set(["utils", "helpers", "misc", "tmp"]);

const TYPOGRAPHIC_POINTS = [0x2018, 0x2019, 0x201c, 0x201d, 0x2013, 0x2014, 0xfe0f];

const TYPOGRAPHIC_CLASS = TYPOGRAPHIC_POINTS.map((point) => `\\u{${point.toString(16)}}`).join("");

export const TYPOGRAPHIC = new RegExp(
  `[${TYPOGRAPHIC_CLASS}]|[\\u{1F000}-\\u{1FAFF}]|[\\u{2600}-\\u{27BF}]|[\\u{2B00}-\\u{2BFF}]`,
  "u",
);

export const TEST_CALLERS = new Set(["it", "test", "describe", "bench", "suite"]);

export function calleeName(node: ESTree.CallExpression): string | undefined {
  if (node.callee.type === "Identifier") {
    return node.callee.name;
  }

  if (node.callee.type === "MemberExpression" && node.callee.property.type === "Identifier") {
    return node.callee.property.name;
  }

  return undefined;
}

export function rootCalleeName(node: ESTree.CallExpression): string | undefined {
  let current: ESTree.Expression | ESTree.Super = node.callee;

  while (current.type === "MemberExpression") {
    current = current.object;
  }

  return current.type === "Identifier" ? current.name : undefined;
}

export function memberChain(node: ESTree.CallExpression): string[] {
  const parts: string[] = [];
  let current: ESTree.Expression | ESTree.Super = node.callee;

  while (current.type === "MemberExpression") {
    if (current.property.type === "Identifier") {
      parts.unshift(current.property.name);
    }

    current = current.object;
  }

  if (current.type === "Identifier") {
    parts.unshift(current.name);
  }

  return parts;
}

export function declaredIdentifier(node: ESTree.Node): DeclaredIdentifier | undefined {
  const candidate =
    "id" in node && node.id !== null && node.id !== undefined ? node.id : "key" in node ? node.key : undefined;

  return candidate !== undefined && candidate.type === "Identifier" ? candidate : undefined;
}

export function isNode(value: unknown): value is ESTree.Node {
  return typeof value === "object" && value !== null && "type" in value && typeof value.type === "string";
}

interface DeclaredIdentifier extends Ranged {
  type: "Identifier";
  name: string;
}
