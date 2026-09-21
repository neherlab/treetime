import { noVagueIdentifiersRule } from "../rules/no-vague-identifiers.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

tester.run("treetime/no-vague-identifiers", noVagueIdentifiersRule, {
  valid: ["const alignmentColumns = 3", "function parseNewick() {}", "class SequenceStore {}"],
  invalid: [
    { code: "const utils = {}", errors: [{ messageId: "vague", data: { name: "utils" } }] },
    { code: "function helpers() {}", errors: [{ messageId: "vague", data: { name: "helpers" } }] },
    { code: "import misc from './misc'", errors: [{ messageId: "vague", data: { name: "misc" } }] },
  ],
});
