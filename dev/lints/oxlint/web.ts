import { definePlugin } from "@oxlint/plugins";

import { useThemedCnRule } from "./rules/use-themed-cn.ts";

export default definePlugin({
  meta: { name: "web" },
  rules: {
    "use-themed-cn": useThemedCnRule,
  },
});
