import { useClassNameHelperRule } from "../rules/use-class-name-helper.ts"
import { ruleTester } from "./rule-tester.ts"

const tester = ruleTester("ts")

tester.run("treetime/use-class-name-helper", useClassNameHelperRule, {
  valid: ["import { cn } from '@/ui/cn'", "import { useState } from 'react'"],
  invalid: [
    { code: "import clsx from 'clsx'", errors: [{ messageId: "useCn", data: { source: "clsx" } }] },
    {
      code: "import { twMerge } from 'tailwind-merge'",
      errors: [{ messageId: "useCn", data: { source: "tailwind-merge" } }],
    },
  ],
})
