import { defineRule } from "@oxlint/plugins";

const DISABLE_DIRECTIVE = /^\s*oxlint-disable(?:-next-line|-line)?\b/;

const HAS_REASON = /\s--\s+\S/;

export const requireSuppressionReasonRule = defineRule({
  meta: {
    type: "suggestion",
    docs: {
      description: "Require a ` -- reason` on every oxlint-disable directive.",
    },
    messages: {
      missingReason:
        "State why this suppression is justified after ` -- `. An unexplained disable hides an unaddressed problem.",
    },
  },
  createOnce(context) {
    return {
      Program() {
        for (const comment of context.sourceCode.getAllComments()) {
          if (DISABLE_DIRECTIVE.test(comment.value) && !HAS_REASON.test(comment.value)) {
            context.report({ loc: comment.loc, messageId: "missingReason" });
          }
        }
      },
    };
  },
});
