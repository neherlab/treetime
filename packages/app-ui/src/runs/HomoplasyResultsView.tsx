import type {
  GapFill,
  HomoplasyResults,
  HomoplasySite,
  HomoplasyStatistics,
  RunRecord,
  RunResults,
} from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";
import TriangleAlert from "~icons/lucide/triangle-alert";

import { nucleotideColorBy, nucleotidePosition } from "../auspice/genotype";
import { Panel, SummaryStrip } from "../components/Panel";
import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
import { Button } from "../ui/button";
import { Empty, EmptyDescription } from "../ui/empty";
import { GenomeSitesChart } from "./GenomeSitesChart";
import { drmText, homoplasySummary, initialHomoplasyColorBy, pressedMutations, siteAt } from "./homoplasy";
import { HomoplasyTables } from "./HomoplasyTables";
import { OutputFiles } from "./OutputFiles";
import { RecurrentTable } from "./RecurrentTable";
import { MultiplicityChart, SiteHitsChart } from "./SiteHitsChart";
import { MissingTree, TreeView, type TreeData } from "./TreeView";
import type { TreeLink } from "./TreeWorkspace";

export function HomoplasyResultsView({
  record,
  results,
  data,
  tree,
}: {
  record: RunRecord;
  results: RunResults;
  data: HomoplasyResults;
  tree: TreeData | undefined;
}) {
  const statistics = data.statistics;
  const colorBy = useMemo(() => initialHomoplasyColorBy(statistics), [statistics]);

  const summary = useMemo(
    () => (statistics === undefined ? undefined : homoplasySummary(record, statistics)),
    [record, statistics],
  );

  const gapFill = record.command === "homoplasy" ? record.config.gap_fill : undefined;

  const aside = useCallback(
    (link: TreeLink) => (statistics === undefined ? undefined : <HomoplasyAside statistics={statistics} link={link} />),
    [statistics],
  );

  const below = useCallback(
    (link: TreeLink) =>
      statistics === undefined ? undefined : <HomoplasyPanels statistics={statistics} gapFill={gapFill} link={link} />,
    [gapFill, statistics],
  );

  return (
    <div className="grid gap-4">
      {summary !== undefined && <SummaryStrip entries={summary} />}
      {statistics === undefined && (
        <Alert>
          <TriangleAlert aria-hidden />
          <AlertTitle>No homoplasy statistics</AlertTitle>
          <AlertDescription>
            This run wrote no readable statistics file. The output files are listed below.
          </AlertDescription>
        </Alert>
      )}
      {tree === undefined ? (
        <MissingTree />
      ) : (
        <TreeView data={tree} colorBy={colorBy} entropy={false} aside={aside} below={below} />
      )}
      <OutputFiles record={record} citation={results.citation} />
    </div>
  );
}

function HomoplasyAside({ statistics, link }: { statistics: HomoplasyStatistics; link: TreeLink }) {
  const position = nucleotidePosition(link.colorBy);
  const colorBySite = useColorBySite(link.setColorBy);
  const site = siteAt(statistics.sites, position);
  const pressed = useMemo(() => pressedMutations(statistics.recurrent, position), [position, statistics.recurrent]);

  return (
    <>
      <Panel title="Recurrent mutations" hint="Select a mutation to color the tree by the base at its position">
        {statistics.recurrent.length === 0 ? (
          <Empty className="py-6">
            <EmptyDescription>No substitution occurs on more than one branch.</EmptyDescription>
          </Empty>
        ) : (
          <RecurrentTable
            label="Recurrent mutations"
            rows={statistics.recurrent}
            pressed={pressed}
            drmAnnotated={statistics.drm_annotated}
            onColor={colorBySite}
          />
        )}
      </Panel>
      {site !== undefined && <SelectedPosition site={site} link={link} />}
    </>
  );
}

function SelectedPosition({ site, link }: { site: HomoplasySite; link: TreeLink }) {
  return (
    <Panel
      title={`Position ${site.display_position}, ${site.branches} branches`}
      hint="Letters other than A, C, G, T in the legend are ambiguous base calls"
      actions={
        link.zoomed ? (
          <Button type="button" variant="outline" size="xs" onClick={link.reset}>
            Show whole tree
          </Button>
        ) : undefined
      }
    >
      <div className="grid max-h-[28rem] gap-3 overflow-auto px-4 py-3">
        {site.substitutions.map((substitution) => (
          <section key={substitution.mutation} aria-label={substitution.mutation} className="grid gap-1.5">
            <h3 className="text-sm">
              <span className="font-mono font-bold">{substitution.mutation}</span> on {substitution.branches}{" "}
              {substitution.branches === 1 ? "branch" : "branches"}
              {substitution.drm !== undefined && (
                <span className="text-muted-foreground"> ({drmText(substitution.drm)})</span>
              )}
            </h3>
            <ul className="flex flex-wrap gap-1">
              {substitution.branch_names.map((name) => (
                <li key={name}>
                  <BranchButton name={name} onSelect={link.select} />
                </li>
              ))}
            </ul>
          </section>
        ))}
      </div>
    </Panel>
  );
}

function BranchButton({ name, onSelect }: { name: string; onSelect: (name: string) => void }) {
  const select = useCallback(() => onSelect(name), [name, onSelect]);

  return (
    <Button type="button" variant="outline" size="xs" onClick={select}>
      {name}
    </Button>
  );
}

function HomoplasyPanels({
  statistics,
  gapFill,
  link,
}: {
  statistics: HomoplasyStatistics;
  gapFill: GapFill | undefined;
  link: TreeLink;
}) {
  const position = nucleotidePosition(link.colorBy);
  const colorBySite = useColorBySite(link.setColorBy);

  return (
    <>
      <Panel
        figure
        title="Sites hit more than once along the genome"
        hint="Substitutions between A, C, G and T; ambiguous base calls are not counted"
      >
        {statistics.sites.length === 0 ? (
          <Empty className="py-6">
            <EmptyDescription>No site has substitutions on more than one branch.</EmptyDescription>
          </Empty>
        ) : (
          <GenomeSitesChart
            sites={statistics.sites}
            genomeLength={statistics.genome_length}
            zeroBased={statistics.zero_based}
            selected={position}
            onSelect={colorBySite}
          />
        )}
      </Panel>
      <div className="grid gap-4 @4xl:grid-cols-[minmax(0,2fr)_minmax(0,1fr)]">
        <Panel
          figure
          title="Substitutions per site"
          hint="Bars: observed sites. Line: expected under a Poisson distribution with the same mean"
        >
          <SiteHitsChart rows={statistics.site_hits} />
        </Panel>
        <Panel
          figure
          title="Branches per mutation"
          hint="Distinct substitutions by the number of branches they occur on"
        >
          <MultiplicityChart rows={statistics.multiplicities} />
        </Panel>
      </div>
      <HomoplasyTables
        statistics={statistics}
        gapFill={gapFill}
        link={link}
        position={position}
        onColor={colorBySite}
      />
    </>
  );
}

function useColorBySite(setColorBy: (colorBy: string) => void): (position: number) => void {
  return useCallback(
    (position: number) => {
      const key = nucleotideColorBy(position);

      if (key !== undefined) {
        setColorBy(key);
      }
    },
    [setColorBy],
  );
}
