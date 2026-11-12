import type { CheckLevel, RunCheck, RunSummary } from "@neherlab/app-contracts";
import { runsGet } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@neherlab/design/ui/collapsible";
import { Kbd } from "@neherlab/design/ui/kbd";
import { Spinner } from "@neherlab/design/ui/spinner";
import { formatForDisplay } from "@tanstack/react-hotkeys";
import { Link } from "@tanstack/react-router";
import { use, useCallback } from "react";
import { useFormContext, useFormState } from "react-hook-form";
import CircleAlert from "~icons/lucide/circle-alert";
import CircleCheck from "~icons/lucide/circle-check";
import Copy from "~icons/lucide/copy";
import Info from "~icons/lucide/info";
import OctagonX from "~icons/lucide/octagon-x";
import Play from "~icons/lucide/play";

import { useApi } from "../api/hooks";
import { RUN_HOTKEY } from "../hotkeys";
import { commandSettings } from "../settings/catalog";
import { SEED_KEY } from "../settings/rerun";
import { useDraftStore } from "../store/draft";
import { toFormValue, type FormConfig } from "./formValues";
import { WorkspaceContext } from "./workspaceContext";

const LEVEL_ICON: Record<CheckLevel, React.ReactNode> = {
  block: <OctagonX aria-label="Blocks the run" className="text-destructive size-4" />,
  warn: <CircleAlert aria-label="Warning" className="text-warning size-4" />,
  advice: <Info aria-label="Advice" className="text-muted-foreground size-4" />,
};

export function RunFooter({
  checks,
  duplicate,
  verb,
  onRun,
}: {
  checks: readonly RunCheck[] | undefined;
  duplicate: RunSummary | undefined;
  verb: string;
  onRun: () => void;
}) {
  const { isSubmitting, isValid, errors } = useFormState<FormConfig>();
  const formError = errors.root?.message;
  const checking = checks === undefined;
  const blocking = checks?.some((check) => check.level === "block") ?? false;
  const findings = checks ?? [];
  const first = findings.find((check) => check.level === "block") ?? findings.at(0);

  const reason = checking
    ? "Checking the settings"
    : blocking
      ? "Resolve the blocking checks first"
      : isValid
        ? `Run with ${formatForDisplay(RUN_HOTKEY)}`
        : "Correct the settings marked as invalid first";

  return (
    <div className="grid gap-2">
      <Collapsible className="grid gap-1.5 text-sm">
        <div className="flex min-w-0 items-start gap-1.5">
          {formError === undefined ? (
            <CheckSummary checking={checking} first={first} />
          ) : (
            <>
              <span className="mt-0.5">{LEVEL_ICON.block}</span>
              <span className="min-w-0 flex-1">{formError}</span>
            </>
          )}
          {findings.length > 0 && (
            <CollapsibleTrigger className="text-primary ml-auto shrink-0 text-xs underline-offset-4 hover:underline">
              {findings.length === 1 ? "Details" : `All ${findings.length} checks`}
            </CollapsibleTrigger>
          )}
        </div>
        <CollapsibleContent>
          <ul className="grid max-h-56 gap-1.5 overflow-y-auto overscroll-contain">
            {findings.map((check) => (
              <CheckItem key={check.id} check={check} locate />
            ))}
          </ul>
        </CollapsibleContent>
      </Collapsible>
      {duplicate !== undefined && <DuplicateRun duplicate={duplicate} />}
      <Button
        type="button"
        size="lg"
        className="w-full"
        disabled={checking || blocking || !isValid || isSubmitting}
        onClick={onRun}
      >
        {isSubmitting ? <Spinner /> : <Play aria-hidden />}
        {duplicate === undefined ? verb : `${verb} again`}
        <Kbd className="bg-primary-foreground/15 text-primary-foreground group-disabled/button:bg-foreground/10 group-disabled/button:text-muted-foreground ml-auto">
          {formatForDisplay(RUN_HOTKEY)}
        </Kbd>
      </Button>
      {(checking || blocking || !isValid) && <p className="text-muted-foreground text-center text-xs">{reason}</p>}
    </div>
  );
}

function CheckSummary({ checking, first }: { checking: boolean; first: RunCheck | undefined }) {
  if (checking) {
    return (
      <>
        <Spinner className="mt-0.5 size-4" />
        <span className="text-muted-foreground">Checking the settings</span>
      </>
    );
  }

  if (first === undefined) {
    return (
      <>
        <CircleCheck aria-hidden className="text-success mt-0.5 size-4 shrink-0" />
        <span>Ready to run.</span>
      </>
    );
  }

  return (
    <>
      <span className="mt-0.5">{LEVEL_ICON[first.level]}</span>
      <span className="line-clamp-2 min-w-0 flex-1">{first.text}</span>
    </>
  );
}

export function CheckItem({ check, locate = false }: { check: RunCheck; locate?: boolean }) {
  const { setValue } = useFormContext<FormConfig>();
  const command = use(WorkspaceContext)?.command;
  const fix = check.fix ?? null;
  const [first] = check.settings;

  const target =
    locate && command !== undefined && first !== undefined
      ? commandSettings(command).settings.find((spec) => spec.key === first)
      : undefined;

  const apply = useCallback(() => {
    for (const setting of fix?.patch ?? []) {
      setValue(setting.key, toFormValue(setting.value), {
        shouldDirty: true,
        shouldValidate: true,
      });
    }
  }, [fix, setValue]);

  const show = useCallback(() => {
    if (target !== undefined) {
      useDraftStore.getState().focusSetting(target.key);
    }
  }, [target]);

  return (
    <li className="grid grid-cols-[1.125rem_1fr] items-start gap-1.5">
      <span className="mt-0.5">{LEVEL_ICON[check.level]}</span>
      <span>{check.text}</span>
      {(fix !== null || target !== undefined) && (
        <div className="col-start-2 flex flex-wrap items-center gap-1.5">
          {fix !== null && (
            <Button type="button" variant="outline" size="sm" onClick={apply}>
              {fix.label}
            </Button>
          )}
          {target !== undefined && (
            <Button type="button" variant="ghost" size="sm" onClick={show}>
              Show {target.label}
            </Button>
          )}
        </div>
      )}
    </li>
  );
}

function DuplicateRun({ duplicate }: { duplicate: RunSummary }) {
  const { data: record } = useApi((context) => runsGet({ ...context, path: { id: duplicate.id } }));
  const seed = record?.id === duplicate.id ? record.seed : undefined;
  const chosen = record?.changed_settings.includes(SEED_KEY) ?? false;

  return (
    <p className="text-muted-foreground flex items-start gap-1.5 text-xs">
      <Copy aria-hidden className="mt-px size-3.5 shrink-0" />
      <span>
        Run &quot;{duplicate.title}&quot; has the same inputs and settings.{" "}
        <Link
          to="/runs/$id/results"
          params={{ id: duplicate.id }}
          className="text-primary font-bold underline-offset-4 hover:underline"
        >
          Open it
        </Link>
        .
        {seed !== undefined && (
          <>
            {" "}
            The run has a random step with seed <span className="font-mono">{seed}</span>
            {chosen ? ", which these settings also set, so a new run repeats it." : "; a new run draws another seed."}
          </>
        )}
      </span>
    </p>
  );
}
