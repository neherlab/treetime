import { defineRule } from "@oxlint/plugins";
import type { ESTree } from "@oxlint/plugins";
import * as z from "zod";

const MIN_MEMBERS = 2;

const zEnumNode = z.object({ enum: z.array(z.string()).min(MIN_MEMBERS) });

const zTaggedNode = z.object({ discriminator: z.object({ mapping: z.record(z.string(), z.string()) }) });

const zConstantsNode = z.object({ oneOf: z.array(z.object({ const: z.string() })).min(MIN_MEMBERS) });

const zOptions = z.object({ enums: z.array(z.array(z.string())) });

export const noContractEnumCopyRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description:
        "Disallow string literal sets that copy the values of a generated enum; use the generated type or value list instead.",
    },
    messages: {
      enumCopy:
        "These literals ({{values}}) copy values of a generated enum. Use the generated type or value list, so a new value cannot be missed.",
    },
    schema: [
      {
        type: "object",
        properties: { enums: { type: "array", items: { type: "array", items: { type: "string" } } } },
        additionalProperties: false,
      },
    ],
  },
  createOnce(context) {
    function enums(): ReadonlyArray<ReadonlySet<string>> {
      const options = zOptions.safeParse(context.options[0]);

      return options.success ? options.data.enums.map((values) => new Set(values)) : [];
    }

    function check(node: ESTree.Node, values: readonly string[]): void {
      const distinct = new Set(values);

      if (distinct.size >= MIN_MEMBERS && enums().some((known) => [...distinct].every((value) => known.has(value)))) {
        context.report({ node, messageId: "enumCopy", data: { values: [...distinct].join(", ") } });
      }
    }

    return {
      TSUnionType(node) {
        const values = node.types.map(stringLiteralType);

        if (values.every((value) => value !== undefined)) {
          check(node, values);
        }
      },
      TSAsExpression(node) {
        const isConst =
          node.typeAnnotation.type === "TSTypeReference" &&
          node.typeAnnotation.typeName.type === "Identifier" &&
          node.typeAnnotation.typeName.name === "const";

        const values = node.expression.type === "ArrayExpression" ? node.expression.elements.map(stringLiteral) : [];

        if (isConst && values.length > 0 && values.every((value) => value !== undefined)) {
          check(node, values);
        }
      },
      LogicalExpression(node) {
        if (node.parent.type === "LogicalExpression" && node.parent.operator === node.operator) {
          return;
        }

        const comparisons = flatten(node, node.operator).map((operand) => comparedLiteral(operand, context.sourceCode));

        if (comparisons.some((comparison) => comparison === undefined)) {
          return;
        }

        const operands = new Set(comparisons.map((comparison) => comparison?.operand));

        if (operands.size === 1) {
          check(
            node,
            comparisons.flatMap((comparison) => (comparison === undefined ? [] : [comparison.value])),
          );
        }
      },
    };
  },
});

export function openApiEnums(document: string): string[][] {
  const found: string[][] = [];

  JSON.parse(document, (_key, value: unknown) => {
    const values = enumValues(value);

    if (values !== undefined) {
      found.push(values);
    }

    return value;
  });

  return found;
}

function enumValues(value: unknown): string[] | undefined {
  const direct = zEnumNode.safeParse(value);

  if (direct.success) {
    return direct.data.enum;
  }

  const tagged = zTaggedNode.safeParse(value);

  if (tagged.success) {
    return Object.keys(tagged.data.discriminator.mapping);
  }

  const constants = zConstantsNode.safeParse(value);

  return constants.success ? constants.data.oneOf.map((branch) => branch.const) : undefined;
}

function stringLiteralType(type: ESTree.TSType): string | undefined {
  return type.type === "TSLiteralType" && type.literal.type === "Literal" && typeof type.literal.value === "string"
    ? type.literal.value
    : undefined;
}

function stringLiteral(element: ESTree.ArrayExpression["elements"][number]): string | undefined {
  return element?.type === "Literal" && typeof element.value === "string" ? element.value : undefined;
}

function flatten(node: ESTree.Expression, operator: string): ESTree.Expression[] {
  return node.type === "LogicalExpression" && node.operator === operator
    ? [...flatten(node.left, operator), ...flatten(node.right, operator)]
    : [node];
}

function comparedLiteral(
  node: ESTree.Expression,
  sourceCode: { getText(node: ESTree.Node): string },
): { operand: string; value: string } | undefined {
  if (node.type !== "BinaryExpression" || (node.operator !== "===" && node.operator !== "!==")) {
    return undefined;
  }

  const right = stringLiteral(node.right.type === "Literal" ? node.right : null);
  const left = stringLiteral(node.left.type === "Literal" ? node.left : null);

  if (right !== undefined && node.left.type !== "Literal") {
    return { operand: sourceCode.getText(node.left), value: right };
  }

  if (left !== undefined && node.right.type !== "Literal") {
    return { operand: sourceCode.getText(node.right), value: left };
  }

  return undefined;
}
