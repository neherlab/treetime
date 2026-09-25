import type { AppCommand, RunSummaryResult } from "@neherlab/app-contracts";
import { Link } from "@tanstack/react-router";
import { CircleAlert, CircleCheck, Info, OctagonX } from "lucide-react";
import { useCallback } from "react";
import { useFormContext, useFormState } from "react-hook-form";

import { COMMAND_SETTINGS } from "../settings/catalog";
import { hasBlockingCheck, type Check, type CheckLevel } from "../settings/checks";
import type { JsonObject } from "../settings/json";
import { autoTitle } from "../settings/titles";
import { useDraftStore } from "../store/draft";
import { Button } from "../ui";
import { toFormValue, type FormConfig } from "./formValues";

const LEVEL_ICON: Record<CheckLevel, React.ReactNode> = {
  block: <OctagonX size={15} aria-label="Blocks the run" className="text-signal-danger" />,
  warn: <CircleAlert size={15} aria-label="Warning" className="text-signal-warn" />,
  advice: <Info size={15} aria-label="Advice" className="text-ink-faint" />,
};

export function ChecksPanel({
  command,
  config,
  checks,
  duplicate,
  verb,
}: {
  command: AppCommand;
  config: JsonObject;
  checks: readonly Check[];
  duplicate: RunSummaryResult | undefined;
  verb: string;
}) {
  const { isSubmitting, isValid } = useFormState<FormConfig>();
  const title = useDraftStore((state) => state.title);
  const update = useDraftStore((state) => state.update);
  const blocking = hasBlockingCheck(checks);

  const onTitle = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => update({ title: event.target.value }),
    [update],
  );

  const reason = blocking
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
          {checks.map((check) => (
            <CheckItem key={check.id} check={check} />
          ))}
          {checks.length === 0 && (
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
        <label className="text-ink-muted grid gap-1">
          Run name
          <input
            type="text"
            value={title}
            onChange={onTitle}
            placeholder={autoTitle(command, COMMAND_SETTINGS[command].specs, config)}
            className="border-line-strong bg-surface-1 text-ink rounded-md border px-2.5 py-1.5"
          />
        </label>
        <Button type="submit" className="w-full" disabled={blocking || !isValid || isSubmitting} title={reason}>
          {duplicate === undefined ? verb : `${verb} again`}
        </Button>
        {(blocking || !isValid) && <p className="text-ink-faint text-center text-xs">{reason}</p>}
      </div>
    </div>
  );
}

function CheckItem({ check }: { check: Check }) {
  const { setValue } = useFormContext<FormConfig>();
  const fix = check.fix;

  const apply = useCallback(() => {
    if (fix !== null) {
      setValue(fix.path.join("."), toFormValue(fix.value), { shouldDirty: true, shouldValidate: true });
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
