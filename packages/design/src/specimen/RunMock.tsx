import { CartesianGrid, ReferenceLine, Scatter, ScatterChart, XAxis, YAxis } from "recharts";
import Check from "~icons/lucide/check";
import Download from "~icons/lucide/download";
import Plus from "~icons/lucide/plus";
import RefreshCw from "~icons/lucide/refresh-cw";
import TriangleAlert from "~icons/lucide/triangle-alert";

import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
import { Badge } from "../ui/badge";
import { Button } from "../ui/button";
import { Card, CardHeader, CardTitle } from "../ui/card";
import { ChartContainer } from "../ui/chart";
import { DescriptionDetails, DescriptionList, DescriptionTerm } from "../ui/description-list";
import { Stat } from "../ui/stat";
import { Wordmark } from "../ui/wordmark";
import { TableDemo } from "./displays";

const CLOCK_RATE = 0.003_54;

const ROOT_DATE = 2009.4;

const TIPS = Array.from({ length: 40 }, (_, index) => {
  const date = 2011.5 + index * 0.22;

  return { date, divergence: CLOCK_RATE * (date - ROOT_DATE) + 0.0012 * Math.sin(index * 1.7) };
});

const OUTLIERS = [
  { date: 2014.2, divergence: 0.0355 },
  { date: 2016.9, divergence: 0.0105 },
];

const FIT_SEGMENT = [
  { x: ROOT_DATE, y: 0 },
  { x: 2020.4, y: CLOCK_RATE * (2020.4 - ROOT_DATE) },
] as const;

const X_DOMAIN = [2009, 2021];

const Y_DOMAIN = [0, 0.045];

const CHART_MARGIN = { top: 8, right: 8, bottom: 0, left: 0 };

const RUNS = [
  { name: "flu/h3n2/200 timetree", state: "Finished in 14 s", current: true },
  { name: "ebola/362 ancestral", state: "Iteration 2 of 3", current: false },
  { name: "zika/86 clock", state: "Finished with 1 warning", current: false },
  { name: "mpox/clade-ii/500 timetree", state: "The tree file cannot be read", current: false },
] as const;

export function RunMock() {
  return (
    <div className="bg-background overflow-hidden rounded-xl border">
      <header className="bg-card flex h-12 items-center gap-3 border-b px-4">
        <Wordmark name="treetime" />
        <span className="flex-1" />
        <Button size="sm">
          <Plus aria-hidden />
          New run
        </Button>
      </header>
      <div className="grid @3xl:grid-cols-[14rem_minmax(0,1fr)]">
        <nav aria-label="Runs" className="bg-card hidden flex-col gap-1 border-r p-2 @3xl:flex">
          {RUNS.map((run) => (
            <a
              key={run.name}
              href="#run"
              aria-current={run.current ? "page" : undefined}
              className="hover:bg-muted aria-[current=page]:bg-selected flex flex-col rounded-md px-2 py-1.5"
            >
              <span className="text-sm font-medium">{run.name}</span>
              <span className="text-muted-foreground text-xs">{run.state}</span>
            </a>
          ))}
        </nav>
        <RunPage />
      </div>
    </div>
  );
}

function RunPage() {
  return (
    <div className="flex min-w-0 flex-col gap-4 p-4">
      <div className="flex flex-wrap items-center gap-3">
        <h3 className="text-lg font-semibold">flu/h3n2/200 timetree</h3>
        <Badge variant="success">
          <Check aria-hidden />
          Finished
        </Badge>
        <span className="flex-1" />
        <Button variant="outline">
          <Download aria-hidden />
          Download
        </Button>
        <Button>
          <RefreshCw aria-hidden />
          Run again
        </Button>
      </div>
      <Card size="sm" className="grid grid-cols-2 gap-4 p-4 @3xl:grid-cols-4">
        <Stat label="Clock rate" value="3.54" uncertainty="± 0.21 × 10⁻³" />
        <Stat label="Root date" value="2009.4" uncertainty="2008.9 to 2009.8" />
        <Stat label="R² of the regression" value="0.92" />
        <Stat label="Tips" value="200" uncertainty="3 without date" />
      </Card>
      <Alert variant="warning">
        <TriangleAlert aria-hidden />
        <AlertTitle>3 sequences have no date</AlertTitle>
        <AlertDescription>
          They stay in the tree, and TreeTime infers their dates; the clock fit leaves them out.
        </AlertDescription>
      </Alert>
      <div className="grid gap-4 @4xl:grid-cols-[3fr_2fr]">
        <RegressionCard />
        <Card size="sm">
          <CardHeader>
            <CardTitle>Settings</CardTitle>
          </CardHeader>
          <DescriptionList className="px-4">
            <DescriptionTerm>Tree</DescriptionTerm>
            <DescriptionDetails className="font-mono">tree.nwk</DescriptionDetails>
            <DescriptionTerm>Alignment</DescriptionTerm>
            <DescriptionDetails className="font-mono">aln.fasta.xz</DescriptionDetails>
            <DescriptionTerm>Model</DescriptionTerm>
            <DescriptionDetails>GTR, inferred</DescriptionDetails>
            <DescriptionTerm>Coalescent</DescriptionTerm>
            <DescriptionDetails>Constant, optimized</DescriptionDetails>
            <DescriptionTerm>Seed</DescriptionTerm>
            <DescriptionDetails className="font-mono">1749203</DescriptionDetails>
          </DescriptionList>
        </Card>
      </div>
      <TableDemo />
    </div>
  );
}

function RegressionCard() {
  return (
    <Card size="sm">
      <CardHeader>
        <CardTitle>Root-to-tip regression</CardTitle>
      </CardHeader>
      <ChartContainer className="aspect-[2/1] px-2" aria-label="Root-to-tip divergence against sampling date">
        <ScatterChart margin={CHART_MARGIN}>
          <CartesianGrid stroke="var(--chart-grid)" vertical={false} />
          <XAxis type="number" dataKey="date" domain={X_DOMAIN} stroke="var(--chart-axis)" tickLine={false} />
          <YAxis type="number" dataKey="divergence" domain={Y_DOMAIN} stroke="var(--chart-axis)" tickLine={false} />
          <ReferenceLine segment={FIT_SEGMENT} stroke="var(--chart-fit)" strokeWidth={1.5} />
          <Scatter data={TIPS} fill="var(--chart-series)" isAnimationActive={false} />
          <Scatter
            data={OUTLIERS}
            fill="none"
            stroke="var(--chart-fault)"
            strokeWidth={1.5}
            isAnimationActive={false}
          />
        </ScatterChart>
      </ChartContainer>
    </Card>
  );
}
