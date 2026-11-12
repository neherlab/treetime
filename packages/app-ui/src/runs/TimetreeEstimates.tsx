import type { CoalescentPrior, RunRecord, SparseConfig, TimetreeEstimates } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@neherlab/design/ui/collapsible";
import { DescriptionDetails, DescriptionList, DescriptionTerm } from "@neherlab/design/ui/description-list";
import { Spinner } from "@neherlab/design/ui/spinner";
import { Stat } from "@neherlab/design/ui/stat";
import type { ReactNode } from "react";
import ChevronRight from "~icons/lucide/chevron-right";
import TriangleAlert from "~icons/lucide/triangle-alert";

import { formatDuration, formatLevel, yearDateText } from "../format";
import { formatSignificant, RATE_UNIT, valueWithSpread } from "../results/estimates";
import { fromJsonFloat, nonFiniteLabel } from "../results/numbers";
import { warningsOf } from "./RunWarnings";
import { useStartVariant } from "./useCreateVariant";

const CORRELATION_DIGITS = 3;

const LIKELIHOOD_DECIMALS = 1;

const WITH_CONFIDENCE: SparseConfig = { confidence: true };

export function TimetreeEstimatesPanel({ record, estimates }: { record: RunRecord; estimates: TimetreeEstimates }) {
  return (
    <section aria-label="Estimates" className="bg-card grid gap-3 rounded-lg border px-4 py-3">
      <div className="flex flex-wrap items-start gap-x-10 gap-y-3">
        <RootDate record={record} estimates={estimates} />
        <ClockRate estimates={estimates} />
        <dl className="grid grid-cols-[max-content_minmax(0,1fr)] gap-x-4 gap-y-1 self-center text-sm">
          <dt className="text-muted-foreground">Temporal signal</dt>
          <dd className="tabular-nums">{temporalSignalText(estimates)}</dd>
          <dt className="text-muted-foreground">Samples used</dt>
          <dd>
            <SamplesUsed record={record} estimates={estimates} />
          </dd>
        </dl>
      </div>
      <RunDetails record={record} estimates={estimates} />
    </section>
  );
}

function RootDate({ record, estimates }: { record: RunRecord; estimates: TimetreeEstimates }) {
  const date = estimates.root_date;
  const interval = estimates.root_interval;
  const edge = estimates.root_near_interval_edge;

  return (
    <Stat
      label="Root date"
      value={date === undefined ? "not dated" : yearDateText(date)}
      detail={
        interval === undefined ? (
          <NoInterval record={record} />
        ) : (
          <>
            <span className="tabular-nums">
              {formatLevel(interval.level)} interval: {yearDateText(interval.lower)} to {yearDateText(interval.upper)}
            </span>
            {edge && <Caution>The estimate lies at the edge of its interval</Caution>}
          </>
        )
      }
    />
  );
}

function NoInterval({ record }: { record: RunRecord }) {
  const { start, busy } = useStartVariant(record, WITH_CONFIDENCE);

  return (
    <span className="flex flex-wrap items-center gap-x-2">
      No interval computed
      <Button type="button" variant="link" size="xs" className="h-auto p-0" onClick={start} disabled={busy}>
        {busy && <Spinner />}
        Compute the interval
      </Button>
    </span>
  );
}

function ClockRate({ estimates }: { estimates: TimetreeEstimates }) {
  const rate = estimates.clock_rate;

  if (rate === undefined) {
    return <Stat label="Clock rate" value="not written" />;
  }

  const { value, spread } = valueWithSpread(rate, estimates.clock_rate_std);

  return (
    <Stat
      label="Clock rate"
      value={value}
      uncertainty={spread}
      detail={estimates.clock_rate_fixed ? `${RATE_UNIT}, fixed` : RATE_UNIT}
    />
  );
}

function SamplesUsed({ record, estimates }: { record: RunRecord; estimates: TimetreeEstimates }) {
  const used = estimates.samples - estimates.excluded_samples;
  const warnings = warningsOf(record.warnings, "samples");

  return (
    <span className="grid gap-0.5">
      <span className="tabular-nums">
        {used} of {estimates.samples}
      </span>
      {estimates.excluded_samples > 0 && (
        <span className="text-muted-foreground text-xs">
          {estimates.excluded_samples} without a usable date or clock outliers
        </span>
      )}
      {warnings.map((warning) => (
        <Caution key={warning.message}>{warning.message}</Caution>
      ))}
    </span>
  );
}

function RunDetails({ record, estimates }: { record: RunRecord; estimates: TimetreeEstimates }) {
  const relax = estimates.relaxed_clock;

  return (
    <Collapsible className="border-t pt-2">
      <CollapsibleTrigger className="text-muted-foreground hover:text-foreground group flex items-center gap-1 text-xs">
        <ChevronRight aria-hidden className="size-3.5 transition-transform group-data-panel-open:rotate-90" />
        Run details
      </CollapsibleTrigger>
      <CollapsibleContent>
        <DescriptionList className="pt-2 text-xs">
          <DescriptionTerm>Coalescent prior</DescriptionTerm>
          <DescriptionDetails>{coalescentPriorText(estimates.coalescent_prior)}</DescriptionDetails>
          <DescriptionTerm>Clock model</DescriptionTerm>
          <DescriptionDetails>
            {relax === undefined ? "Strict clock" : `Relaxed clock: slack ${relax.slack}, coupling ${relax.coupling}`}
          </DescriptionDetails>
          <DescriptionTerm>Log likelihood</DescriptionTerm>
          <DescriptionDetails className="tabular-nums">{likelihoodText(estimates)}</DescriptionDetails>
          <DescriptionTerm>Iterations</DescriptionTerm>
          <DescriptionDetails className="tabular-nums">{estimates.iterations}</DescriptionDetails>
          {record.duration_seconds !== undefined && (
            <>
              <DescriptionTerm>Duration</DescriptionTerm>
              <DescriptionDetails className="tabular-nums">
                {formatDuration(record.duration_seconds)}
              </DescriptionDetails>
            </>
          )}
        </DescriptionList>
      </CollapsibleContent>
    </Collapsible>
  );
}

function Caution({ children }: { children: ReactNode }) {
  return (
    <span className="text-warning flex items-start gap-1 text-xs font-bold">
      <TriangleAlert aria-label="Caution" className="mt-px size-3 shrink-0" />
      <span>{children}</span>
    </span>
  );
}

function temporalSignalText(estimates: TimetreeEstimates): string {
  if (estimates.clock_rate_fixed) {
    return "Fixed clock rate";
  }

  const parts = [
    estimates.r === undefined ? undefined : `r = ${formatSignificant(estimates.r, CORRELATION_DIGITS)}`,
    estimates.r_squared === undefined
      ? undefined
      : `R² = ${formatSignificant(estimates.r_squared, CORRELATION_DIGITS)}`,
  ].filter((part) => part !== undefined);

  return parts.length === 0 ? "not computed" : parts.join(", ");
}

function likelihoodText(estimates: TimetreeEstimates): string {
  if (estimates.log_likelihood === undefined) {
    return "not written";
  }

  const value = fromJsonFloat(estimates.log_likelihood);

  return Number.isFinite(value) ? value.toFixed(LIKELIHOOD_DECIMALS) : nonFiniteLabel(value);
}

export function coalescentPriorText(prior: CoalescentPrior): string {
  if (prior.kind === "fixed") {
    return prior.unit === "substitutions-per-site"
      ? `Constant size, Tc = ${prior.tc} substitutions per site`
      : `Constant size, Tc = ${prior.tc}`;
  }

  if (prior.kind === "optimized") {
    return "Constant size, optimized Tc";
  }

  if (prior.kind === "skyline") {
    return `Skyline, ${prior.points} points, stiffness ${prior.stiffness}`;
  }

  return "None";
}
