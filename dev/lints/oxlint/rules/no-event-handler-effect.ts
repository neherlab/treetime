import { defineRule } from "@oxlint/plugins";

import { hasNonEmptyDeps, makeEffectState, SUBSCRIPTION_CALLEES } from "./effects.ts";

export const noEventHandlerEffectRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow an effect that only forwards to a single event handler call.",
    },
    messages: {
      handler: "This effect only forwards to an event handler. Call the handler from the event instead of an effect.",
    },
  },
  createOnce(context) {
    const { state, visitors } = makeEffectState();

    return {
      ...visitors,
      "Program:exit"() {
        for (const effect of state.effects) {
          const body = effect.callback.body;

          if (
            body === null ||
            !hasNonEmptyDeps(effect.deps) ||
            body.type !== "BlockStatement" ||
            body.body.length !== 1
          ) {
            continue;
          }

          const statement = body.body[0];

          if (
            statement?.type === "ExpressionStatement" &&
            statement.expression.type === "CallExpression" &&
            statement.expression.callee.type === "Identifier" &&
            !state.setterNames.has(statement.expression.callee.name) &&
            !SUBSCRIPTION_CALLEES.has(statement.expression.callee.name)
          ) {
            context.report({ node: effect.node, messageId: "handler" });
          }
        }
      },
    };
  },
});
