import { noTypographicCharactersRule } from "../rules/no-typographic-characters.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts");

const curlyQuote = String.fromCodePoint(0x201c);

const emDash = String.fromCodePoint(0x2014);

tester.run("treetime/no-typographic-characters", noTypographicCharactersRule, {
  valid: ['const label = "straight"', "const dash = 'a-b'"],
  invalid: [
    { code: `const label = "${curlyQuote}quoted"`, errors: [{ messageId: "typographic" }] },
    { code: `const range = "1${emDash}2"`, errors: [{ messageId: "typographic" }] },
  ],
});
