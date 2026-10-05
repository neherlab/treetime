import type { CheckLevel, RunCheck, RunSummary } from "@neherlab/app-contracts";
import { formatForDisplay } from "@tanstack/react-hotkeys";
import { Link } from "@tanstack/react-router";
import { useCallback } from "react";
import { useFormContext, useFormState } from "react-hook-form";
import CircleAlert from "~icons/lucide/circle-alert";
import CircleCheck from "~icons/lucide/circle-check";
import Copy from "~icons/lucide/copy";
import Info from "~icons/lucide/info";
import OctagonX from "~icons/lucide/octagon-x";
import Play from "~icons/lucide/play";

import { Panel } from "../components/Panel";
import { RUN_HOTKEY } from "../hotkeys";
import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
import { Button } from "../ui/button";
import { Kbd } from "../ui/kbd";
import { Spinner } from "../ui/spinner";
import { toFormValue, type FormConfig } from "./formValues";

const LEVEL_ICON: Record<CheckLevel, React.ReactNode> = {
  block: <OctagonX aria-label="Blocks the run" className="text-destructive size-4" />,
  warn: <CircleAlert aria-label="Warning" className="text-warning size-4" />,
  advice: <Info aria-label="Advice" className="text-muted-foreground size-4" />,
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
        ? `Run with ${formatForDisplay(RUN_HOTKEY)}`
        : "Correct the settings marked as invalid first";

  return (
    <Panel title="Checks">
      <div className="grid gap-3 p-3.5">
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
              <CircleCheck aria-hidden className="text-success mt-0.5 size-4" />
              <span>Ready to run.</span>
            </li>
          )}
        </ul>
        {duplicate !== undefined && (
          <Alert className="border-warning/40 bg-warning/10">
            <Copy aria-hidden />
            <AlertTitle>Run &quot;{duplicate.title}&quot; has the same inputs and settings</AlertTitle>
            <AlertDescription>
              <span>
                <Link
                  to="/runs/$id/results"
                  params={{ id: duplicate.id }}
                  className="text-primary font-bold underline-offset-4 hover:underline"
                >
                  Open it
                </Link>{" "}
                or run again to check that it reproduces.
              </span>
            </AlertDescription>
          </Alert>
        )}
        <Button
          type="submit"
          size="lg"
          className="w-full"
          disabled={checking || blocking || !isValid || isSubmitting}
          title={reason}
        >
          {isSubmitting ? <Spinner /> : <Play aria-hidden />}
          {duplicate === undefined ? verb : `${verb} again`}
          <Kbd className="bg-primary-foreground/15 text-primary-foreground ml-auto">{formatForDisplay(RUN_HOTKEY)}</Kbd>
        </Button>
        {(checking || blocking || !isValid) && <p className="text-muted-foreground text-center text-xs">{reason}</p>}
      </div>
    </Panel>
  );
}

export function CheckItem({ check }: { check: RunCheck }) {
  const { setValue } = useFormContext<FormConfig>();
  const fix = check.fix ?? null;

  const apply = useCallback(() => {
    for (const setting of fix?.patch ?? []) {
      setValue(setting.path.join("."), toFormValue(setting.value), {
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
