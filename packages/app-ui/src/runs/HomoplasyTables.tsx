import type { AmbiguousSite, GapFill, HomoplasyStatistics, TaxonResult } from "@neherlab/app-contracts";
import { memo, useCallback, useMemo } from "react";
import ChevronRight from "~icons/lucide/chevron-right";

import { DataTable, dataColumns } from "../components/DataTable";
import { Panel } from "../components/Panel";
import { Button } from "../ui/button";
import { Card } from "../ui/card";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "../ui/collapsible";
import { Empty, EmptyDescription } from "../ui/empty";
import { gapFillNote, pressedMutations } from "./homoplasy";
import { PositionButton, RecurrentTable } from "./RecurrentTable";
import type { TreeLink } from "./TreeWorkspace";

const taxonColumn = dataColumns<TaxonResult>();

const ambiguousColumn = dataColumns<AmbiguousSite>();

const TAXON_NUMERIC = new Set(["homoplasic", "drm", "ambiguous", "indels"]);

const AMBIGUOUS_NUMERIC = new Set(["branches"]);

const HOMOPLASIC_SORT = [{ id: "homoplasic", desc: true }];

const BRANCHES_SORT = [{ id: "branches", desc: true }];

export function HomoplasyTables({
  statistics,
  gapFill,
  link,
  position,
  onColor,
}: {
  statistics: HomoplasyStatistics;
  gapFill: GapFill | undefined;
  link: TreeLink;
  position: number | undefined;
  onColor: (position: number) => void;
}) {
  const pressedIndels = useMemo(
    () => pressedMutations(statistics.recurrent_indels, position),
    [position, statistics.recurrent_indels],
  );

  return (
    <>
      <Panel
        title="Samples with homoplasic substitutions"
        hint="Substitutions on the terminal branch at sites hit more than once; select a sample to find it in the tree"
      >
        <TaxaTable taxa={statistics.taxa} drmAnnotated={statistics.drm_annotated} onSelect={link.select} />
      </Panel>
      <Panel
        title="Recurrent insertions and deletions"
        hint="Select one to color the tree by the base at its first column"
      >
        {statistics.recurrent_indels.length === 0 ? (
          <Empty className="py-6">
            <EmptyDescription>No insertion or deletion occurs on more than one branch.</EmptyDescription>
          </Empty>
        ) : (
          <RecurrentTable
            label="Recurrent insertions and deletions"
            rows={statistics.recurrent_indels}
            pressed={pressedIndels}
            drmAnnotated={false}
            onColor={onColor}
          />
        )}
      </Panel>
      <AmbiguousPanel statistics={statistics} gapFill={gapFill} position={position} onColor={onColor} />
    </>
  );
}

const TaxaTable = memo(function TaxaTable({
  taxa,
  drmAnnotated,
  onSelect,
}: {
  taxa: readonly TaxonResult[];
  drmAnnotated: boolean;
  onSelect: (name: string) => void;
}) {
  const columns = useMemo(() => taxonColumns(drmAnnotated, onSelect), [drmAnnotated, onSelect]);

  return (
    <DataTable
      label="Samples with homoplasic substitutions"
      columns={columns}
      rows={taxa}
      rowId={taxonKey}
      initialSorting={HOMOPLASIC_SORT}
      numeric={TAXON_NUMERIC}
    />
  );
});

function SampleButton({ name, onSelect }: { name: string; onSelect: (name: string) => void }) {
  const select = useCallback(() => onSelect(name), [name, onSelect]);

  return (
    <Button type="button" variant="ghost" size="xs" onClick={select}>
      {name}
    </Button>
  );
}

function AmbiguousPanel({
  statistics,
  gapFill,
  position,
  onColor,
}: {
  statistics: HomoplasyStatistics;
  gapFill: GapFill | undefined;
  position: number | undefined;
  onColor: (position: number) => void;
}) {
  const note = gapFillNote(gapFill);
  const shown = statistics.ambiguous_sites.length;

  const columns = useMemo(() => ambiguousColumns(position, onColor), [onColor, position]);

  return (
    <Collapsible>
      <Card size="sm" className="min-w-0 gap-0 py-0">
        <CollapsibleTrigger className="group/trigger hover:bg-muted/50 flex w-full items-center gap-2 px-4 py-2.5 text-left font-bold">
          <ChevronRight
            aria-hidden
            className="text-muted-foreground size-4 transition-transform group-data-panel-open/trigger:rotate-90"
          />
          Ambiguous characters ({statistics.ambiguous_changes} changes)
        </CollapsibleTrigger>
        <CollapsibleContent className="grid gap-2 border-t px-4 py-3">
          {note !== undefined && <p className="text-muted-foreground text-sm">{note}</p>}
          <p className="text-muted-foreground text-xs">
            Top {shown} of {statistics.ambiguous_site_count} sites
          </p>
          <DataTable
            label="Sites with ambiguous changes"
            columns={columns}
            rows={statistics.ambiguous_sites}
            rowId={ambiguousKey}
            initialSorting={BRANCHES_SORT}
            numeric={AMBIGUOUS_NUMERIC}
          />
        </CollapsibleContent>
      </Card>
    </Collapsible>
  );
}

function taxonColumns(drmAnnotated: boolean, onSelect: (name: string) => void) {
  return taxonColumn.columns([
    taxonColumn.accessor((row) => row.name, {
      id: "sample",
      header: "Sample",
      cell: ({ getValue }) => <SampleButton name={getValue()} onSelect={onSelect} />,
    }),
    taxonColumn.accessor((row) => row.homoplasic_mutations.length, { id: "homoplasic", header: "Homoplasic" }),
    ...(drmAnnotated ? [taxonColumn.accessor((row) => row.drm_mutations ?? 0, { id: "drm", header: "DRM" })] : []),
    taxonColumn.accessor((row) => row.ambiguous_changes, { id: "ambiguous", header: "Ambiguous changes" }),
    taxonColumn.accessor((row) => row.recurrent_indels, { id: "indels", header: "Recurrent indels" }),
    taxonColumn.accessor((row) => row.homoplasic_mutations.join(" "), {
      id: "mutations",
      header: "Mutations",
      cell: ({ getValue }) => <span className="font-mono whitespace-normal">{getValue()}</span>,
    }),
  ]);
}

function ambiguousColumns(position: number | undefined, onColor: (position: number) => void) {
  return ambiguousColumn.columns([
    ambiguousColumn.accessor((row) => row.display_position, {
      id: "position",
      header: "Position",
      cell: ({ row }) => (
        <PositionButton
          position={row.original.position}
          text={String(row.original.display_position)}
          label={`Color the tree by the base at position ${row.original.display_position}`}
          pressed={row.original.position === position}
          onColor={onColor}
        />
      ),
    }),
    ambiguousColumn.accessor((row) => row.branches, { id: "branches", header: "Branches" }),
  ]);
}

function taxonKey(row: TaxonResult): string {
  return row.name;
}

function ambiguousKey(row: AmbiguousSite): string {
  return String(row.position);
}
