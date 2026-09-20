import { noStateInEffectInitializerRule } from "../rules/no-state-in-effect-initializer.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("tsx")

tester.run("web/no-state-in-effect-initializer", noStateInEffectInitializerRule, {
  valid: [
    "const [count, setCount] = useState(0); useEffect(() => { document.title = String(count) }, [count])",
  ],
  invalid: [
    {
      code: "const [count, setCount] = useState(0); useEffect(() => { setCount(compute()) }, [])",
      errors: [{ messageId: "initializer" }],
    },
  ],
})
