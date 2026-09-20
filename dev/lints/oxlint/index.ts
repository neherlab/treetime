import { definePlugin } from "@oxlint/plugins"

import { noAssertionInLoopRule } from "./rules/no-assertion-in-loop.ts"
import { noAsyncArrayPredicateRule } from "./rules/no-async-array-predicate.ts"
import { noDisabledTestsRule } from "./rules/no-disabled-tests.ts"
import { noExecShellStringRule } from "./rules/no-exec-shell-string.ts"
import { noFakeSuccessRule } from "./rules/no-fake-success.ts"
import { noFocusedTestsRule } from "./rules/no-focused-tests.ts"
import { noLeakyMocksRule } from "./rules/no-leaky-mocks.ts"
import { noModuleLevelMutableRule } from "./rules/no-module-level-mutable.ts"
import { noTautologicalAssertionRule } from "./rules/no-tautological-assertion.ts"
import { noTestResourceAccessRule } from "./rules/no-test-resource-access.ts"
import { noTypographicCharactersRule } from "./rules/no-typographic-characters.ts"
import { noUnboundedSuppressionRule } from "./rules/no-unbounded-suppression.ts"
import { noUppercaseTestTitleRule } from "./rules/no-uppercase-test-title.ts"
import { noVagueIdentifiersRule } from "./rules/no-vague-identifiers.ts"
import { noVersionedNamesRule } from "./rules/no-versioned-names.ts"
import { preferStrictEqualRule } from "./rules/prefer-strict-equal.ts"
import { preferTestOverItRule } from "./rules/prefer-test-over-it.ts"
import { requireIoTimeoutRule } from "./rules/require-io-timeout.ts"
import { requireSuppressionReasonRule } from "./rules/require-suppression-reason.ts"
import { useClassNameHelperRule } from "./rules/use-class-name-helper.ts"

export default definePlugin({
  meta: { name: "treetime" },
  rules: {
    "no-assertion-in-loop": noAssertionInLoopRule,
    "no-async-array-predicate": noAsyncArrayPredicateRule,
    "no-disabled-tests": noDisabledTestsRule,
    "no-exec-shell-string": noExecShellStringRule,
    "no-fake-success": noFakeSuccessRule,
    "no-focused-tests": noFocusedTestsRule,
    "no-leaky-mocks": noLeakyMocksRule,
    "no-module-level-mutable": noModuleLevelMutableRule,
    "no-tautological-assertion": noTautologicalAssertionRule,
    "no-test-resource-access": noTestResourceAccessRule,
    "no-typographic-characters": noTypographicCharactersRule,
    "no-unbounded-suppression": noUnboundedSuppressionRule,
    "no-uppercase-test-title": noUppercaseTestTitleRule,
    "no-vague-identifiers": noVagueIdentifiersRule,
    "no-versioned-names": noVersionedNamesRule,
    "prefer-strict-equal": preferStrictEqualRule,
    "prefer-test-over-it": preferTestOverItRule,
    "require-io-timeout": requireIoTimeoutRule,
    "require-suppression-reason": requireSuppressionReasonRule,
    "use-class-name-helper": useClassNameHelperRule,
  },
})
