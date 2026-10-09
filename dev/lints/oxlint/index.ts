import { definePlugin } from "@oxlint/plugins";

import { callersBeforeCalleesRule } from "./rules/callers-before-callees.ts";
import { noAssertionInLoopRule } from "./rules/no-assertion-in-loop.ts";
import { noContractAliasRule } from "./rules/no-contract-alias.ts";
import { noContractEnumCopyRule } from "./rules/no-contract-enum-copy.ts";
import { noContractParseRule } from "./rules/no-contract-parse.ts";
import { noDeclareGlobalRule } from "./rules/no-declare-global.ts";
import { noExecShellStringRule } from "./rules/no-exec-shell-string.ts";
import { noFakeSuccessRule } from "./rules/no-fake-success.ts";
import { noLeakyMocksRule } from "./rules/no-leaky-mocks.ts";
import { noModuleLevelMutableRule } from "./rules/no-module-level-mutable.ts";
import { noSideEffectsInGettersRule } from "./rules/no-side-effects-in-getters.ts";
import { noTautologicalAssertionRule } from "./rules/no-tautological-assertion.ts";
import { noTestResourceAccessRule } from "./rules/no-test-resource-access.ts";
import { noTypographicCharactersRule } from "./rules/no-typographic-characters.ts";
import { noUnboundedSuppressionRule } from "./rules/no-unbounded-suppression.ts";
import { noVagueIdentifiersRule } from "./rules/no-vague-identifiers.ts";
import { noVersionedNamesRule } from "./rules/no-versioned-names.ts";
import { requireIoTimeoutRule } from "./rules/require-io-timeout.ts";
import { requireSuppressionReasonRule } from "./rules/require-suppression-reason.ts";
import { useThemedCnRule } from "./rules/use-themed-cn.ts";

export default definePlugin({
  meta: { name: "custom" },
  rules: {
    "callers-before-callees": callersBeforeCalleesRule,
    "no-assertion-in-loop": noAssertionInLoopRule,
    "no-contract-alias": noContractAliasRule,
    "no-contract-enum-copy": noContractEnumCopyRule,
    "no-contract-parse": noContractParseRule,
    "no-declare-global": noDeclareGlobalRule,
    "no-exec-shell-string": noExecShellStringRule,
    "no-fake-success": noFakeSuccessRule,
    "no-leaky-mocks": noLeakyMocksRule,
    "no-module-level-mutable": noModuleLevelMutableRule,
    "no-side-effects-in-getters": noSideEffectsInGettersRule,
    "no-tautological-assertion": noTautologicalAssertionRule,
    "no-test-resource-access": noTestResourceAccessRule,
    "no-typographic-characters": noTypographicCharactersRule,
    "no-unbounded-suppression": noUnboundedSuppressionRule,
    "no-vague-identifiers": noVagueIdentifiersRule,
    "no-versioned-names": noVersionedNamesRule,
    "require-io-timeout": requireIoTimeoutRule,
    "require-suppression-reason": requireSuppressionReasonRule,
    "use-themed-cn": useThemedCnRule,
  },
});
