import { toNestErrors } from "@hookform/resolvers";
import {
  zAncestralConfig,
  zClockConfig,
  zHomoplasyConfig,
  zMugrationConfig,
  zOptimizeConfig,
  zPruneConfig,
  zTimetreeConfig,
  type AppCommand,
} from "@neherlab/app-contracts";
import type { FieldError, Resolver } from "react-hook-form";

import type { FormConfig } from "./formValues";

const ROOT_ERROR = "root";

const CONFIG_SCHEMAS = {
  timetree: zTimetreeConfig,
  optimize: zOptimizeConfig,
  prune: zPruneConfig,
  ancestral: zAncestralConfig,
  homoplasy: zHomoplasyConfig,
  clock: zClockConfig,
  mugration: zMugrationConfig,
} as const;

export function configResolver(command: AppCommand): Resolver<FormConfig> {
  const schema = CONFIG_SCHEMAS[command];

  return async (values, _context, options) => {
    const result = await schema.safeParseAsync(values);

    if (result.success) {
      return { values, errors: {} };
    }

    const errors: Record<string, FieldError> = {};

    for (const issue of result.error.issues) {
      const path = issue.path.length === 0 ? ROOT_ERROR : issue.path.map(String).join(".");

      errors[path] ??= { type: issue.code, message: issue.message };
    }

    return { values: {}, errors: toNestErrors(errors, options) };
  };
}
