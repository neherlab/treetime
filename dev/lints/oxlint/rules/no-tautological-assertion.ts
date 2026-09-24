import { defineRule } from "@oxlint/plugins";
import type { ESTree, SourceCode, Variable } from "@oxlint/plugins";

import { bindingImport, bindingImportName, bindingVariable } from "./binding.ts";

const UNKNOWN = { kind: "unknown" } as const;

const EXPECT_EQUALITY = new Set(["toBe", "toEqual", "toStrictEqual"]);

const EXPECT_MODULES = new Set(["vitest", "@jest/globals"]);

const ASSERT_MODULES = new Set(["assert", "node:assert", "assert/strict", "node:assert/strict"]);

const ASSERT_EQUALITY = new Map<string, Comparison>([
  ["equal", { mode: "nodeLoose", negated: false }],
  ["notEqual", { mode: "nodeLoose", negated: true }],
  ["strictEqual", { mode: "sameValue", negated: false }],
  ["notStrictEqual", { mode: "sameValue", negated: true }],
  ["deepEqual", { mode: "deepLoose", negated: false }],
  ["notDeepEqual", { mode: "deepLoose", negated: true }],
  ["deepStrictEqual", { mode: "deepStrict", negated: false }],
  ["notDeepStrictEqual", { mode: "deepStrict", negated: true }],
]);

const BINARY_EQUALITY = new Map<string, Comparison>([
  ["==", { mode: "loose", negated: false }],
  ["!=", { mode: "loose", negated: true }],
  ["===", { mode: "strict", negated: false }],
  ["!==", { mode: "strict", negated: true }],
]);

export const noTautologicalAssertionRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow assertions proved to pass from constant values or stable identity.",
    },
    messages: {
      tautology:
        "This assertion is proved to pass without checking a result. Compare against an independent expected value.",
    },
  },
  createOnce(context) {
    return {
      CallExpression(node) {
        const result = expectResult(node, context.sourceCode) ?? assertResult(node, context.sourceCode);

        if (result === true) {
          context.report({ node, messageId: "tautology" });
        }
      },
    };
  },
});

function expectResult(node: ESTree.CallExpression, sourceCode: SourceCode): boolean | undefined {
  const { callee } = node;

  if (callee.type !== "MemberExpression") {
    return undefined;
  }

  const matcher = memberName(callee);
  let subject = callee.object;
  let negated = false;

  if (subject.type === "MemberExpression" && memberName(subject) === "not") {
    negated = true;
    subject = subject.object;
  }

  if (
    subject.type !== "CallExpression" ||
    subject.callee.type !== "Identifier" ||
    assertionBinding(subject.callee, sourceCode)?.kind !== "expect"
  ) {
    return undefined;
  }

  const actual = argumentAt(subject, 0);

  if (actual === undefined || matcher === undefined) {
    return undefined;
  }

  if (matcher === "toBeTruthy" || matcher === "toBeFalsy") {
    const truth = truthValue(expressionFact(actual, sourceCode));

    return applyPolarity(truth, negated !== (matcher === "toBeFalsy"));
  }

  const expected = argumentAt(node, 0);

  if (!EXPECT_EQUALITY.has(matcher) || expected === undefined) {
    return undefined;
  }

  const mode = matcher === "toBe" ? "sameValue" : "deepStrict";

  return applyPolarity(
    compareFacts(expressionFact(actual, sourceCode), expressionFact(expected, sourceCode), mode),
    negated,
  );
}

function assertResult(node: ESTree.CallExpression, sourceCode: SourceCode): boolean | undefined {
  const { callee } = node;
  const actual = argumentAt(node, 0);

  if (actual === undefined) {
    return undefined;
  }

  const identifier = callee.type === "MemberExpression" ? callee.object : callee;

  if (identifier.type !== "Identifier") {
    return undefined;
  }

  const binding = assertionBinding(identifier, sourceCode);

  if (binding?.kind !== "assert" || (callee.type === "Identifier" && !binding.callable)) {
    return undefined;
  }

  const method = callee.type === "MemberExpression" ? memberName(callee) : binding.method;

  if (method === "ok" || (callee.type === "Identifier" && method === undefined)) {
    return truthValue(expressionFact(actual, sourceCode));
  }

  const comparison = method === undefined ? undefined : ASSERT_EQUALITY.get(method);
  const expected = argumentAt(node, 1);

  if (comparison === undefined || expected === undefined) {
    return undefined;
  }

  const left = expressionFact(actual, sourceCode);
  const right = expressionFact(expected, sourceCode);
  const result = compareFacts(left, right, comparison.mode);

  if (comparison.mode === "nodeLoose" || comparison.mode === "deepLoose") {
    const strictMode = comparison.mode === "nodeLoose" ? "sameValue" : "deepStrict";
    const strictResult = compareFacts(left, right, strictMode);

    if (binding.strict === true) {
      return applyPolarity(strictResult, comparison.negated);
    }

    if (binding.strict === undefined && result !== strictResult) {
      return undefined;
    }
  }

  return applyPolarity(result, comparison.negated);
}

function assertionBinding(
  identifier: ESTree.IdentifierReference,
  sourceCode: SourceCode,
): AssertionBinding | undefined {
  const binding = bindingImport(identifier, sourceCode);

  if (binding === undefined) {
    const variable = bindingVariable(identifier, sourceCode);

    if (variable !== undefined && variable.defs.length > 0) {
      return undefined;
    }

    if (identifier.name === "expect") {
      return { kind: "expect" };
    }

    return identifier.name === "assert"
      ? { kind: "assert", method: undefined, strict: undefined, callable: true }
      : undefined;
  }

  const name = bindingImportName(binding);

  if (EXPECT_MODULES.has(binding.source) && name === "expect") {
    return { kind: "expect" };
  }

  if (!ASSERT_MODULES.has(binding.source)) {
    return undefined;
  }

  return {
    kind: "assert",
    method: name === "strict" || name === "default" ? undefined : name,
    strict: binding.source.endsWith("/strict") || name === "strict",
    callable: binding.node.type !== "ImportNamespaceSpecifier",
  };
}

function argumentAt(node: ESTree.CallExpression, index: number): ESTree.Expression | undefined {
  const argument = node.arguments[index];

  return argument?.type === "SpreadElement" ? undefined : argument;
}

function memberName(node: ESTree.MemberExpression): string | undefined {
  return !node.computed && node.property.type === "Identifier" ? node.property.name : undefined;
}

function expressionFact(
  expression: ESTree.Expression,
  sourceCode: SourceCode,
  seen: ReadonlySet<Variable> = new Set(),
): Fact {
  if (expression.type === "Literal") {
    return "regex" in expression
      ? { kind: "reference", identity: expression }
      : { kind: "primitive", value: expression.value };
  }

  if (expression.type === "Identifier") {
    return identifierFact(expression, sourceCode, seen);
  }

  if (expression.type === "UnaryExpression") {
    return unaryFact(expression, sourceCode, seen);
  }

  if (expression.type === "BinaryExpression" && expression.left.type !== "PrivateIdentifier") {
    const left = expressionFact(expression.left, sourceCode, seen);
    const right = expressionFact(expression.right, sourceCode, seen);
    const comparison = BINARY_EQUALITY.get(expression.operator);

    const result =
      comparison === undefined
        ? relationalResult(left, right, expression.operator)
        : applyPolarity(compareFacts(left, right, comparison.mode), comparison.negated);

    return result === undefined ? UNKNOWN : { kind: "primitive", value: result };
  }

  if (expression.type === "TemplateLiteral" && expression.expressions.length === 0) {
    const value = expression.quasis[0]?.value.cooked;

    return value === undefined || value === null ? UNKNOWN : { kind: "primitive", value };
  }

  if (expression.type === "ObjectExpression" || expression.type === "ArrayExpression") {
    return { kind: "reference", identity: expression };
  }

  if (
    expression.type === "TSAsExpression" ||
    expression.type === "TSSatisfiesExpression" ||
    expression.type === "TSNonNullExpression" ||
    expression.type === "ParenthesizedExpression"
  ) {
    return expressionFact(expression.expression, sourceCode, seen);
  }

  return UNKNOWN;
}

function identifierFact(
  identifier: ESTree.IdentifierReference,
  sourceCode: SourceCode,
  seen: ReadonlySet<Variable>,
): Fact {
  const variable = bindingVariable(identifier, sourceCode);

  if (variable === undefined || variable.defs.length === 0) {
    if (sourceCode.isGlobalReference(identifier)) {
      if (identifier.name === "NaN") {
        return { kind: "primitive", value: NaN };
      }

      if (identifier.name === "undefined") {
        return { kind: "primitive", value: undefined };
      }

      if (identifier.name === "Infinity") {
        return { kind: "primitive", value: Infinity };
      }
    }

    return UNKNOWN;
  }

  if (seen.has(variable) || variable.references.some((reference) => reference.isWrite() && !reference.init)) {
    return UNKNOWN;
  }

  const definition = variable.defs[0];

  if (definition?.type === "ImportBinding" && definition.node.type !== "ImportNamespaceSpecifier") {
    return UNKNOWN;
  }

  if (
    (definition?.type === "Variable" || definition?.type === "ClassName") &&
    definition.node.end >= identifier.start
  ) {
    return UNKNOWN;
  }

  if (
    definition?.type === "Variable" &&
    definition.node.type === "VariableDeclarator" &&
    definition.node.id.type === "Identifier" &&
    definition.node.init !== null &&
    definition.parent?.type === "VariableDeclaration" &&
    definition.parent.kind !== "var" &&
    definition.node.end < identifier.start
  ) {
    const value = expressionFact(definition.node.init, sourceCode, new Set([...seen, variable]));

    if (value.kind !== "unknown") {
      return value;
    }
  }

  return { kind: "binding", variable };
}

function unaryFact(expression: ESTree.UnaryExpression, sourceCode: SourceCode, seen: ReadonlySet<Variable>): Fact {
  const argument = expressionFact(expression.argument, sourceCode, seen);

  if (expression.operator === "!") {
    const truth = truthValue(argument);

    return truth === undefined ? UNKNOWN : { kind: "primitive", value: !truth };
  }

  if (argument.kind !== "primitive") {
    return UNKNOWN;
  }

  if (expression.operator === "-" && (typeof argument.value === "number" || typeof argument.value === "bigint")) {
    return { kind: "primitive", value: -argument.value };
  }

  if (expression.operator === "+" && typeof argument.value === "number") {
    return argument;
  }

  if (expression.operator === "void") {
    return { kind: "primitive", value: undefined };
  }

  return UNKNOWN;
}

function compareFacts(left: Fact, right: Fact, mode: EqualityMode): boolean | undefined {
  if (left.kind === "primitive" && right.kind === "primitive") {
    return comparePrimitives(left.value, right.value, mode);
  }

  if (left.kind === "reference" && right.kind === "reference") {
    if (left.identity === right.identity) {
      return true;
    }

    return mode === "deepLoose" || mode === "deepStrict" ? undefined : false;
  }

  if (left.kind === "binding" && right.kind === "binding" && left.variable === right.variable) {
    return mode === "strict" || mode === "loose" ? undefined : true;
  }

  return undefined;
}

function comparePrimitives(left: Primitive, right: Primitive, mode: EqualityMode): boolean {
  if (mode === "strict") {
    return left === right;
  }

  if (mode === "sameValue" || mode === "deepStrict") {
    return Object.is(left, right);
  }

  if (mode !== "loose" && Number.isNaN(left) && Number.isNaN(right)) {
    return true;
  }

  // oxlint-disable-next-line eslint/eqeqeq -- models the loose equality of assert.equal and assert.notEqual
  return left == right;
}

function relationalResult(left: Fact, right: Fact, operator: string): boolean | undefined {
  if (left.kind !== "primitive" || right.kind !== "primitive") {
    return undefined;
  }

  if (
    (typeof left.value !== "number" || typeof right.value !== "number") &&
    (typeof left.value !== "string" || typeof right.value !== "string")
  ) {
    return undefined;
  }

  switch (operator) {
    case "<":
      return left.value < right.value;
    case "<=":
      return left.value <= right.value;
    case ">":
      return left.value > right.value;
    case ">=":
      return left.value >= right.value;
    default:
      return undefined;
  }
}

function truthValue(fact: Fact): boolean | undefined {
  if (fact.kind === "primitive") {
    return Boolean(fact.value);
  }

  return fact.kind === "reference" ? true : undefined;
}

function applyPolarity(result: boolean | undefined, negated: boolean): boolean | undefined {
  return result === undefined ? undefined : result !== negated;
}

interface Comparison {
  mode: EqualityMode;
  negated: boolean;
}

type EqualityMode = "strict" | "loose" | "sameValue" | "nodeLoose" | "deepStrict" | "deepLoose";

type Primitive = string | number | bigint | boolean | null | undefined;

type Fact =
  | { kind: "unknown" }
  | { kind: "primitive"; value: Primitive }
  | { kind: "reference"; identity: ESTree.Node }
  | { kind: "binding"; variable: Variable };

type AssertionBinding =
  | { kind: "expect" }
  | { kind: "assert"; method: string | undefined; strict: boolean | undefined; callable: boolean };
