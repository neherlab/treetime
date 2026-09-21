import { defineRule } from "@oxlint/plugins";

import { effectBodyFacts, hasEmptyDeps, makeEffectState } from "./effects.ts";

export const noStateInEffectInitializerRule = defineRule({
  meta: {
    type: "problem",
    docs: {
      description: "Disallow initializing state inside a mount effect; pass the initial value to useState.",
    },
    messages: {
      initializer: "Do not initialize state inside a mount effect. Pass the initial value to useState instead.",
    },
  },
  createOnce(context) {
    const { state, visitors } = makeEffectState();

    return {
      ...visitors,
      "Program:exit"() {
        for (const effect of state.effects) {
          if (!hasEmptyDeps(effect.deps)) {
            continue;
          }

          const { setterCalls } = effectBodyFacts(effect.callback, state.setterNames);

          if (setterCalls.length > 0) {
            context.report({ node: effect.node, messageId: "initializer" });
          }
        }
      },
    };
  },
});
