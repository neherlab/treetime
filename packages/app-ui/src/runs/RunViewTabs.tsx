import type { RunRecord, RunResults, TimetreeResults } from "@neherlab/app-contracts";
import { Alert, AlertDescription, AlertTitle } from "@neherlab/design/ui/alert";
import { DescriptionDetails, DescriptionList, DescriptionTerm } from "@neherlab/design/ui/description-list";
import { Navigate } from "@tanstack/react-router";
import type { ReactNode } from "react";

import { GENOTYPE_COLORING } from "../auspice/toolbar";
import { LoadingState } from "../components/PageShell";
import { Panel } from "../components/Panel";
import { useLocale } from "../i18n/locale";
import { commandSettings } from "../settings/catalog";
import { ClockView } from "./ClockView";
import { OutputDownload } from "./OutputFiles";
import { useRunOutputs, useRunView } from "./runContext";
import { settingText } from "./settingText";
import { SkylinePlot } from "./SkylinePlot";
import { SubstitutionModels } from "./SubstitutionModels";
import { coalescentPriorText } from "./TimetreeEstimates";
import { TreeView, type TreeData } from "./TreeView";
import { hasRunView, type RunViewKey } from "./views";

const COALESCENT_SETTING = "coalescent";

export function RunClockTab() {
  return (
    <ResultsGate view="clock">
      {(record, results, tree) =>
        results.results.command === "timetree" && (
          <ClockView
            record={record}
            regression={results.results.data.root_to_tip}
            statistics={results.results.data.statistics}
            tree={tree}
          />
        )
      }
    </ResultsGate>
  );
}

export function RunPopulationTab() {
  return (
    <ResultsGate view="population">
      {(record, results) =>
        results.results.command === "timetree" && <PopulationView record={record} data={results.results.data} />
      }
    </ResultsGate>
  );
}

export function RunModelTab() {
  return (
    <ResultsGate view="model">
      {(record, results, tree) => <ModelView record={record} results={results} tree={tree} />}
    </ResultsGate>
  );
}

export function ResultsGate({
  view,
  children,
}: {
  view: RunViewKey;
  children: (record: RunRecord, results: RunResults, tree: TreeData | undefined) => ReactNode;
}) {
  const { record } = useRunView();

  if (!hasRunView(record.command, record.status, view)) {
    return <Navigate to="/runs/$id/results" params={{ id: record.id }} replace />;
  }

  return <ReadyOutputs record={record}>{children}</ReadyOutputs>;
}

function ReadyOutputs({
  record,
  children,
}: {
  record: RunRecord;
  children: (record: RunRecord, results: RunResults, tree: TreeData | undefined) => ReactNode;
}) {
  const outputs = useRunOutputs(record);

  if (outputs.state === "failed") {
    return (
      <Alert variant="destructive">
        <AlertTitle>The outputs of the run cannot be read</AlertTitle>
        <AlertDescription>{outputs.message}</AlertDescription>
      </Alert>
    );
  }

  if (outputs.state === "loading") {
    return <LoadingState text="Reading the outputs" />;
  }

  return children(record, outputs.results, outputs.tree);
}

function PopulationView({ record, data }: { record: RunRecord; data: TimetreeResults }) {
  const locale = useLocale();
  const specs = commandSettings(record.command).settings;
  const group = specs.find((spec) => spec.key === COALESCENT_SETTING)?.group;
  const settings = specs.filter((spec) => spec.group === group && spec.role === "setting");

  return (
    <div className="grid gap-4 @4xl:grid-cols-[minmax(0,2fr)_minmax(0,1fr)]">
      <Panel
        figure
        title="Coalescent time scale"
        hint={
          <>
            Estimate of <i>T</i>
            <sub>c</sub> in years with its interval. The effective population size equals <i>T</i>
            <sub>c</sub> multiplied by the number of generations per year. Shown up to the latest sample.
          </>
        }
        actions={
          <OutputDownload runId={record.id} kind="coalescent">
            Skyline table
          </OutputDownload>
        }
      >
        {data.skyline.length === 0 ? (
          <p className="text-muted-foreground px-2 py-3 text-sm">
            This run estimated no skyline. The coalescent prior of the run is listed beside.
          </p>
        ) : (
          <SkylinePlot segments={data.skyline} until={data.estimates?.latest_sample?.year} />
        )}
      </Panel>
      <Panel title="Coalescent settings">
        <DescriptionList className="px-3.5 py-3">
          {data.estimates !== undefined && (
            <>
              <DescriptionTerm>Prior used</DescriptionTerm>
              <DescriptionDetails>{coalescentPriorText(data.estimates.coalescent_prior)}</DescriptionDetails>
            </>
          )}
          {settings.map((spec) => (
            <div key={spec.key} className="contents">
              <DescriptionTerm>{spec.label}</DescriptionTerm>
              <DescriptionDetails className="font-mono text-xs">
                {settingText(record, spec.key, locale)}
              </DescriptionDetails>
            </div>
          ))}
        </DescriptionList>
      </Panel>
    </div>
  );
}

function ModelView({ record, results, tree }: { record: RunRecord; results: RunResults; tree: TreeData | undefined }) {
  return (
    <div className="grid gap-4">
      {results.models.models.length === 0 ? (
        <Alert>
          <AlertDescription>This run fitted no substitution model.</AlertDescription>
        </Alert>
      ) : (
        <SubstitutionModels record={record} models={results.models} />
      )}
      {tree !== undefined && tree.tree.colorings.some((coloring) => coloring.key === GENOTYPE_COLORING) && (
        <Panel
          title="Entropy"
          hint="Diversity at each position of the alignment; click a bar to color the tree by that position"
        >
          <div className="p-2">
            <TreeView data={tree} colorBy={undefined} entropy />
          </div>
        </Panel>
      )}
    </div>
  );
}
