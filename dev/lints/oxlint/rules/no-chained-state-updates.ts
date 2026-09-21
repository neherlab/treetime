import { defineRule } from "@oxlint/plugins";

import { effectBodyFacts, hasNonEmptyDeps, makeEffectState } from "./effects.ts";

export const noChainedStateUpdatesRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow an effect whose only work is chaining state setters off a dependency change.",
    },
    messages: {
      chained:
        "Do not chain state updates in an effect. Derive the value during render or update it in the event handler.",
    },
  },
  createOnce(context) {
    const { state, visitors } = makeEffectState();

    return {
      ...visitors,
      "Program:exit"() {
        for (const effect of state.effects) {
          if (!hasNonEmptyDeps(effect.deps)) {
            continue;
          }

          const { setterCalls, hasCleanup, statements } = effectBodyFacts(effect.callback, state.setterNames);

          if (hasCleanup || setterCalls.length === 0 || statements.length === 0) {
            continue;
          }

          const onlyStateUpdates = statements.every(
            (statement) =>
              statement.type === "ExpressionStatement" &&
              statement.expression.type === "CallExpression" &&
              statement.expression.callee.type === "Identifier" &&
              state.setterNames.has(statement.expression.callee.name),
          );

          if (onlyStateUpdates) {
            context.report({ node: effect.node, messageId: "chained" });
          }
        }
      },
    };
  },
});
