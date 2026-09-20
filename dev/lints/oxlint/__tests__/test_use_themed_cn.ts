import { useThemedCnRule } from "../rules/use-themed-cn.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("tsx")

tester.run("web/use-themed-cn", useThemedCnRule, {
  valid: ["import { cn } from '@/ui/cn'", "import { useState } from 'react'"],
  invalid: [
    { code: "import clsx from 'clsx'", errors: [{ messageId: "themedCn", data: { source: "clsx" } }] },
    {
      code: "import { cn } from 'cn'",
      errors: [{ messageId: "themedCn", data: { source: "cn" } }],
    },
  ],
})
