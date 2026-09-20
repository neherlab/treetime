import { noFakeSuccessRule } from "../rules/no-fake-success.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("ts", { Promise: "readonly" });
const swallowed = { messageId: "swallowedError" };

tester.run("no_fake_success_binding", noFakeSuccessRule, {
  valid: [
    "Promise.resolve().catch(error => { switch (2) { case 1: let handle; case 2: handle; return [] } })",
    "Promise.resolve().catch(error => { class Holder { [Holder]() {} } return [] })",
    "Promise.resolve().catch(error => { try { pending(log(error)); var pending } finally { return [] } })",
    "Promise.resolve().catch(error => { try { consume(pending, log(error)); if (false) { var pending } } finally { return [] } })",
    "Promise.resolve().catch(error => { try { class Holder { static { consume(pending, log(error)); var pending } } } finally { return [] } })",
    "Promise.resolve().catch(error => { try { switch (1) { case 1: let handle; case 2: handle(log(error)) } } finally { return [] } })",
    "Promise.resolve().catch(error => { try { switch (2) { case 1: let handle; case 2: { let handle; handle(log(error)) } } } finally { return [] } })",
    "Promise.resolve().catch(error => { try { switch (1) { default: let handle; case 2: handle(log(error)) } } finally { return [] } })",
    "Promise.resolve().catch(error => { try { switch (1) { case 1: class Handler {} case 2: new Handler(log(error)) } } finally { return [] } })",
    "Promise.resolve().catch(error => { try { switch (2) { case 1: var pending; case 2: consume(pending, log(error)) } } finally { return [] } })",
    "Promise.resolve().catch(error => { try { let handle; handle(log(error)) } finally { return [] } })",
    "Promise.resolve().catch(error => { try { const handle = () => {}; handle(log(error)) } finally { return [] } })",
    "function run(handle) { return Promise.resolve().catch(error => { try { handle(log(error)) } finally { return [] } }) }",
    "Promise.resolve().catch(error => { try { throw 1 } catch (cause) { try { consume(cause, log(error)) } finally { return [] } } })",
    "Promise.resolve().catch(error => { try { for (let handle; true;) { handle(log(error)); break } } finally { return [] } })",
    "Promise.resolve().catch(error => { switch (1) { case 1: let handle; switch (2) { case 2: handle(log(error)) } } return [] })",
  ],
  invalid: [
    {
      code: "class Holder { [Promise.resolve().catch(error => { Holder; return [] })]() {} }",
      errors: [swallowed],
    },
    {
      code: "Promise.resolve().catch(error => { try { pending(); log(error); var pending } finally { return [] } })",
      errors: [swallowed],
    },
    {
      code: "Promise.resolve().catch(error => { try { pending.method(log(error)); var pending } finally { return [] } })",
      errors: [swallowed],
    },
    {
      code: "Promise.resolve().catch(error => { try { let pending = load(); consume(+pending, log(error)) } finally { return [] } })",
      errors: [swallowed],
    },
    {
      code: "function run(choice) { return Promise.resolve().catch(error => { try { switch (choice) { case 1: let handle; case 2: handle(log(error)); break; default: log(error) } } finally { return [] } }) }",
      errors: [swallowed],
    },
    {
      code: "Promise.resolve().catch(error => { try { switch (2) { case 1: const handle = () => {}; case 2: handle(log(error)) } } finally { return [] } })",
      errors: [swallowed],
    },
    {
      code: "Promise.resolve().catch(error => { try { switch (2) { case 1: class Handler {} case 2: new Handler(log(error)) } } finally { return [] } })",
      errors: [swallowed],
    },
    {
      code: "Promise.resolve().catch(error => { try { switch (2) { default: let handle; case 2: handle(log(error)) } } finally { return [] } })",
      errors: [swallowed],
    },
    {
      code: "Promise.resolve().catch(error => { try { switch (2) { case 1: let handle; case handle(log(error)): break } } finally { return [] } })",
      errors: [swallowed],
    },
    {
      code: "Promise.resolve().catch(error => { try { switch (2) { case 1: let handle; case 2: switch (1) { case 1: handle(log(error)) } } } finally { return [] } })",
      errors: [swallowed],
    },
  ],
});
