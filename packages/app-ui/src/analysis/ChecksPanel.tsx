import type { CheckLevel, RunCheck, RunSummary } from "@neherlab/app-contracts";
import { Link } from "@tanstack/react-router";
import { CircleAlert, CircleCheck, Info, OctagonX } from "lucide-react";
import { useCallback } from "react";
import { useFormContext, useFormState } from "react-hook-form";

import { zJsonValue } from "../settings/json";
import { Button } from "../ui";
import { toFormValue, type FormConfig } from "./formValues";

const LEVEL_ICON: Record<CheckLevel, React.ReactNode> = {
  block: <OctagonX size={15} aria-label="Blocks the run" className="text-signal-danger" />,
  warn: <CircleAlert size={15} aria-label="Warning" className="text-signal-warn" />,
  advice: <Info size={15} aria-label="Advice" className="text-ink-faint" />,
};

export function ChecksPanel({
  checks,
  duplicate,
  verb,
}: {
  checks: readonly RunCheck[] | undefined;
  duplicate: RunSummary | undefined;
  verb: string;
}) {
  const { isSubmitting, isValid, errors } = useFormState<FormConfig>();
  const formError = errors.root?.message;
  const checking = checks === undefined;
  const blocking = checks?.some((check) => check.level === "block") ?? false;

  const reason = checking
    ? "Checking the settings"
    : blocking
      ? "Resolve the blocking checks first"
      : isValid
        ? "Ctrl Enter"
        : "Correct the settings marked as invalid first";

  return (
    <div className="border-line bg-surface-1 rounded-lg border">
      <div className="border-line border-b px-3.5 py-2.5">
        <h3 className="font-bold">Checks</h3>
      </div>
      <div className="grid gap-2.5 px-3.5 py-3">
        <ul className="grid gap-1.5">
          {formError !== undefined && (
            <li className="grid grid-cols-[1.125rem_1fr] items-start gap-1.5">
              <span className="mt-0.5">{LEVEL_ICON.block}</span>
              <span>{formError}</span>
            </li>
          )}
          {(checks ?? []).map((check) => (
            <CheckItem key={check.id} check={check} />
          ))}
          {checks?.length === 0 && (
            <li className="grid grid-cols-[1.125rem_1fr] items-start gap-1.5">
              <CircleCheck size={15} aria-hidden className="text-signal-ok mt-0.5" />
              <span>Ready to run.</span>
            </li>
          )}
        </ul>
        {duplicate !== undefined && (
          <div className="bg-signal-warn-subtle grid gap-1.5 rounded-md px-2.5 py-2">
            <span>Run &quot;{duplicate.title}&quot; has the same inputs and settings.</span>
            <span>
              <Link to="/runs/$id/results" params={{ id: duplicate.id }} className="text-accent font-bold">
                Open it
              </Link>{" "}
              <span className="text-ink-muted">or run again to check that it reproduces</span>
            </span>
          </div>
        )}
        <Button
          type="submit"
          className="w-full"
          disabled={checking || blocking || !isValid || isSubmitting}
          title={reason}
        >
          {duplicate === undefined ? verb : `${verb} again`}
        </Button>
        {(checking || blocking || !isValid) && <p className="text-ink-faint text-center text-xs">{reason}</p>}
      </div>
    </div>
  );
}

export function CheckItem({ check }: { check: RunCheck }) {
  const { setValue } = useFormContext<FormConfig>();
  const fix = check.fix ?? null;

  const apply = useCallback(() => {
    for (const setting of fix?.patch ?? []) {
      setValue(setting.path.join("."), toFormValue(zJsonValue.parse(setting.value)), {
        shouldDirty: true,
        shouldValidate: true,
      });
    }
  }, [fix, setValue]);

  return (
    <li className="grid grid-cols-[1.125rem_1fr] items-start gap-1.5">
      <span className="mt-0.5">{LEVEL_ICON[check.level]}</span>
      <span>{check.text}</span>
      {fix !== null && (
        <Button type="button" variant="outline" size="sm" className="col-start-2 justify-self-start" onClick={apply}>
          {fix.label}
        </Button>
      )}
    </li>
  );
}
