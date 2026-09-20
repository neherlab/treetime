import { noVersionedNamesRule } from "../rules/no-versioned-names.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/no-versioned-names", noVersionedNamesRule, {
  valid: [
    "const parser = {}",
    "function reconstruct() {}",
    "class Cache {}",
  ],
  invalid: [
    { code: "const parser_v2 = {}", errors: [{ messageId: "versioned", data: { name: "parser_v2" } }] },
    { code: "function build_old() {}", errors: [{ messageId: "versioned", data: { name: "build_old" } }] },
    { code: "class Cache_new {}", errors: [{ messageId: "versioned", data: { name: "Cache_new" } }] },
  ],
})
