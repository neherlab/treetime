import { noChainedStateUpdatesRule } from "../rules/no-chained-state-updates.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("tsx");

tester.run("web/no-chained-state-updates", noChainedStateUpdatesRule, {
  valid: ["const [a, setA] = useState(0); useEffect(() => { document.title = String(a) }, [a])"],
  invalid: [
    {
      code: "const [a, setA] = useState(0); const [b, setB] = useState(0); useEffect(() => { setB(a) }, [a])",
      errors: [{ messageId: "chained" }],
    },
  ],
});
