import { useThemedCnRule } from "../rules/use-themed-cn.ts";
import { ruleTester } from "./rule-tester.ts";

const tester = ruleTester("tsx");

tester.run("custom/use-themed-cn", useThemedCnRule, {
  valid: ["import { cn } from '@/ui/cn'", "import { useState } from 'react'", "import clsx from 'clsx'"],
  invalid: [{ code: "import { cn } from 'cn'", errors: [{ messageId: "themedCn" }] }],
});
