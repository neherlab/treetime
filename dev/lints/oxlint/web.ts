import { definePlugin } from "@oxlint/plugins"

import { noChainedStateUpdatesRule } from "./rules/no-chained-state-updates.ts"
import { noEventHandlerEffectRule } from "./rules/no-event-handler-effect.ts"
import { noStateInEffectInitializerRule } from "./rules/no-state-in-effect-initializer.ts"
import { tailwindClassesRule } from "./rules/tailwind-classes.ts"
import { useThemedCnRule } from "./rules/use-themed-cn.ts"

export default definePlugin({
  meta: { name: "web" },
  rules: {
    "no-chained-state-updates": noChainedStateUpdatesRule,
    "no-event-handler-effect": noEventHandlerEffectRule,
    "no-state-in-effect-initializer": noStateInEffectInitializerRule,
    "tailwind-classes": tailwindClassesRule,
    "use-themed-cn": useThemedCnRule,
  },
})
