import { noEventHandlerEffectRule } from "../rules/no-event-handler-effect.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("tsx")

tester.run("web/no-event-handler-effect", noEventHandlerEffectRule, {
  valid: [
    "useEffect(() => { const id = setInterval(tick, 1000); return () => clearInterval(id) }, [tick])",
  ],
  invalid: [
    {
      code: "useEffect(() => { onChange(value) }, [value])",
      errors: [{ messageId: "handler" }],
    },
  ],
})
