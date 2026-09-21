import { defineRule } from "@oxlint/plugins";
import type { ESTree } from "@oxlint/plugins";

import { calleeName, isNode } from "./ast.ts";

export const callersBeforeCalleesRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Order module-level private functions so a caller precedes its callee.",
    },
    messages: {
      outOfOrder:
        "`{{callee}}` is declared before `{{caller}}`, which calls it. Order private functions so a caller precedes its callee.",
    },
  },
  createOnce(context) {
    return {
      Program(program) {
        const position = new Map<string, number>();
        const declaration = new Map<string, ESTree.Node>();
        program.body.forEach((statement, index) => {
          if (statement.type === "FunctionDeclaration" && statement.id !== null) {
            position.set(statement.id.name, index);
            declaration.set(statement.id.name, statement);
          }
        });

        if (position.size < 2) {
          return;
        }

        const edges = new Set<string>();

        for (const statement of program.body) {
          if (statement.type !== "FunctionDeclaration" || statement.id === null) {
            continue;
          }

          const caller = statement.id.name;

          for (const callee of callsWithin(statement, position)) {
            if (callee !== caller) {
              edges.add(`${caller} ${callee}`);
            }
          }
        }

        const reported = new Set<string>();

        for (const edge of edges) {
          const [caller, callee] = edge.split(" ");

          if (caller === undefined || callee === undefined || reported.has(callee)) {
            continue;
          }

          if (edges.has(`${callee} ${caller}`)) {
            continue;
          }

          if ((position.get(caller) ?? 0) > (position.get(callee) ?? 0)) {
            const node = declaration.get(callee);

            if (node !== undefined) {
              reported.add(callee);
              context.report({ node, messageId: "outOfOrder", data: { caller, callee } });
            }
          }
        }
      },
    };
  },
});

function callsWithin(root: ESTree.Node, names: Map<string, number>): Set<string> {
  const found = new Set<string>();
  const stack: ESTree.Node[] = [root];

  while (stack.length > 0) {
    const node = stack.pop();

    if (node === undefined) {
      continue;
    }

    if (node.type === "CallExpression") {
      const name = calleeName(node);

      if (name !== undefined && names.has(name)) {
        found.add(name);
      }
    }

    for (const [key, value] of Object.entries(node)) {
      if (key === "parent") {
        continue;
      }

      if (Array.isArray(value)) {
        for (const item of value) {
          if (isNode(item)) {
            stack.push(item);
          }
        }
      } else if (isNode(value)) {
        stack.push(value);
      }
    }
  }

  return found;
}
