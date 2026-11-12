import { errorMessage, type AppCommand, type InputFacts, type SparseConfig } from "@neherlab/app-contracts";
import { clockPreview } from "@neherlab/app-contracts/client";
import { Alert, AlertDescription, AlertTitle } from "@neherlab/design/ui/alert";
import { Button } from "@neherlab/design/ui/button";
import { Empty, EmptyDescription, EmptyHeader, EmptyTitle } from "@neherlab/design/ui/empty";
import { keepPreviousData } from "@tanstack/react-query";
import { lazy, Suspense, useCallback, useMemo, useState } from "react";
import RefreshCw from "~icons/lucide/refresh-cw";

import { useApi } from "../api/hooks";
import { LoadingState } from "../components/PageShell";
import { Panel } from "../components/Panel";
import { useRootToTip } from "../runs/rootToTip";
import { RootToTipPlot } from "../runs/RootToTipPlot";
import { commandSettings } from "../settings/catalog";
import { inputFactsRequest } from "../settings/inputs";
import { inputFactRows, PREVIEW_ROOT_CAPTIONS } from "./draftPreview";
import { PreviewStatistics } from "./PreviewStatistics";
import { useWorkspace } from "./workspaceContext";

const PreviewTree = lazy(async () => ({ default: (await import("./PreviewTree")).PreviewTree }));

export function DraftMain() {
  const { command, settled, facts } = useWorkspace();
  const hasData = inputFactsRequest(command, settled) !== null;
  const { title } = commandSettings(command);

  if (!hasData) {
    return (
      <Empty className="h-full">
        <EmptyHeader>
          <EmptyTitle>No data yet</EmptyTitle>
          <EmptyDescription>
            Add the input files in the Data section of the options, or open an example dataset.
          </EmptyDescription>
        </EmptyHeader>
      </Empty>
    );
  }

  return (
    <div className="grid content-start gap-4 p-4">
      <header className="grid gap-0.5">
        <h2 className="font-heading text-base font-bold">{title} not run yet</h2>
        <p className="text-muted-foreground text-sm">
          The data as TreeTime reads it. Run the analysis to see its results here.
        </p>
      </header>
      {facts !== undefined && <InputFactsPanel facts={facts} />}
      <ClockPreviewPanel command={command} config={settled} />
    </div>
  );
}

function InputFactsPanel({ facts }: { facts: InputFacts }) {
  const rows = useMemo(() => inputFactRows(facts), [facts]);

  return (
    <Panel title="Input data">
      <dl className="grid grid-cols-[auto_minmax(0,1fr)] gap-x-4 gap-y-1 px-4 py-3 text-sm">
        {rows.map((row) => (
          <div key={row.label} className="contents">
            <dt className="text-muted-foreground">{row.label}</dt>
            <dd className="min-w-0 break-words">{row.value}</dd>
          </div>
        ))}
      </dl>
    </Panel>
  );
}

function ClockPreviewPanel({ command, config }: { command: AppCommand; config: SparseConfig }) {
  const treeKey = JSON.stringify(config["tree"] ?? null);
  const [requestedFor, setRequestedFor] = useState<string | undefined>(undefined);
  const beyondLimit = requestedFor === treeKey;
  const hasTree = config["tree"] !== undefined && config["tree"] !== null;

  const {
    data: preview,
    error,
    isFetching,
  } = useApi(
    (context) =>
      clockPreview({
        ...context,
        body: beyondLimit ? { command, config, beyond_limit: true } : { command, config },
      }),
    { enabled: hasTree, placeholderData: keepPreviousData, staleTime: 30_000, retry: false },
  );

  const update = useCallback(() => setRequestedFor(treeKey), [treeKey]);
  const { points, line } = useRootToTip(preview?.clock?.root_to_tip);
  const run = useMemo(() => ({ id: "preview", title: "Preview", command }), [command]);

  if (!hasTree) {
    return null;
  }

  if (error !== null) {
    return (
      <Alert>
        <AlertTitle>No clock preview</AlertTitle>
        <AlertDescription>{errorMessage(error)}</AlertDescription>
      </Alert>
    );
  }

  if (preview === undefined) {
    return <LoadingState text="Computing the clock preview" />;
  }

  if (preview.tree === undefined) {
    return (
      <Alert>
        <AlertTitle>Preview not computed</AlertTitle>
        <AlertDescription className="grid gap-2">
          <p>
            The tree has {preview.tips} tips. The preview computes trees of up to {preview.tip_limit} tips on every
            change; larger trees on request.
          </p>
          <Button
            type="button"
            variant="outline"
            size="sm"
            className="justify-self-start"
            onClick={update}
            disabled={isFetching}
          >
            <RefreshCw aria-hidden />
            Update preview
          </Button>
        </AlertDescription>
      </Alert>
    );
  }

  return (
    <div className="grid gap-4">
      {preview.clock !== undefined && (
        <Panel figure title="Root-to-tip regression" hint="Dated samples; red points are clock-filter outliers">
          <PreviewStatistics statistics={preview.clock.statistics} />
          <RootToTipPlot points={points} line={line} selected={undefined} inView={undefined} onSelect={undefined} />
        </Panel>
      )}
      <Panel figure title="Tree" hint={PREVIEW_ROOT_CAPTIONS[preview.tree.root]}>
        <Suspense fallback={<LoadingState text="Loading the tree view" />}>
          <PreviewTree document={preview.tree.auspice} run={run} />
        </Suspense>
      </Panel>
    </div>
  );
}
