import { noExecShellStringRule } from "../rules/no-exec-shell-string.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/no-exec-shell-string", noExecShellStringRule, {
  valid: [
    "execFile('git', ['status'])",
    "spawn('ls', ['-la'])",
    "exec(command)",
  ],
  invalid: [
    { code: "exec('rm -rf tmp')", errors: [{ messageId: "execString" }] },
    { code: "execSync(`echo ${value}`)", errors: [{ messageId: "execString" }] },
  ],
})
