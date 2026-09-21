import type { DatasetInfo } from "@neherlab/app-contracts";
import { RotateCcw } from "lucide-react";
import { useCallback, useEffect, useRef } from "react";

import { useDatasets } from "../hooks";
import { useActiveCommand } from "../hooks/useActiveCommand";
import { useAppStore } from "../store/app-store";
import type { CommandName, FileSlotKind } from "../types";
import { FILE_SLOTS } from "../types";
import { Button, cn } from "../ui";
import { FileSlot } from "./FileSlot";

const COMMAND_FILE_REQUIREMENTS: Record<CommandName, { required: FileSlotKind[]; optional: FileSlotKind[] }> = {
  ancestral: { required: ["tree", "alignment"], optional: ["vcfReference"] },
  clock: { required: ["tree", "dates"], optional: ["alignment", "vcfReference"] },
  mugration: { required: ["tree", "states"], optional: ["weights"] },
  optimize: { required: ["tree", "alignment"], optional: [] },
  prune: { required: ["tree"], optional: ["alignment"] },
  timetree: { required: ["tree", "alignment", "dates"], optional: ["vcfReference"] },
};

const QUICK_DATASETS: ReadonlyArray<{ name: string; label: string }> = [
  { name: "flu/h3n2/20", label: "flu20" },
  { name: "flu/h3n2/200", label: "flu200" },
  { name: "ebola/100", label: "ebola100" },
  { name: "sc2/2844", label: "sc2" },
];

const DEFAULT_DATASET = "flu/h3n2/20";

function applyDataset(
  datasetName: string,
  datasets: DatasetInfo[] | undefined,
  setSelectedDataset: (d: string) => void,
  setFile: (kind: FileSlotKind, file: { name: string; size: number } | undefined) => void,
) {
  setSelectedDataset(datasetName);
  const info = datasets?.find((d) => d.name === datasetName);
  const files = info?.files ?? [];

  if (files.includes("tree.nwk")) {
    setFile("tree", { name: `${datasetName}/tree.nwk`, size: 0 });
  }

  if (files.includes("aln.fasta.xz")) {
    setFile("alignment", { name: `${datasetName}/aln.fasta.xz`, size: 0 });
  }

  if (files.includes("metadata.tsv")) {
    setFile("dates", { name: `${datasetName}/metadata.tsv`, size: 0 });
    setFile("states", { name: `${datasetName}/metadata.tsv`, size: 0 });
  }
}

export function FileInputPanel() {
  const activeCommand = useActiveCommand();
  const selectedDataset = useAppStore((s) => s.selectedDataset);
  const setSelectedDataset = useAppStore((s) => s.setSelectedDataset);
  const setFile = useAppStore((s) => s.setFile);
  const resetForm = useAppStore((s) => s.resetForm);
  const reqs = COMMAND_FILE_REQUIREMENTS[activeCommand];
  const { data: datasets } = useDatasets();
  const preloaded = useRef(false);

  useEffect(() => {
    if (datasets && !preloaded.current) {
      preloaded.current = true;
      applyDataset(DEFAULT_DATASET, datasets, setSelectedDataset, setFile);
    }
  }, [datasets, setSelectedDataset, setFile]);

  const handleDatasetChange = useCallback(
    (e: React.ChangeEvent<HTMLSelectElement>) => {
      const dataset = e.target.value;

      if (dataset) {
        applyDataset(dataset, datasets, setSelectedDataset, setFile);
      }
    },
    [setSelectedDataset, setFile, datasets],
  );

  const handleQuickSelect = useCallback(
    (datasetName: string) => {
      applyDataset(datasetName, datasets, setSelectedDataset, setFile);
    },
    [datasets, setSelectedDataset, setFile],
  );

  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <h3 className="text-ink text-sm font-semibold">Input files</h3>
        <Button variant="ghost" size="sm" onClick={resetForm}>
          <RotateCcw size={12} />
          Reset
        </Button>
      </div>

      <div>
        <label htmlFor="dataset-picker" className="text-2xs text-ink-muted mb-1 block">
          Example dataset
        </label>
        <select
          id="dataset-picker"
          value={selectedDataset}
          onChange={handleDatasetChange}
          className="border-line bg-surface-0 text-ink hover:border-line-strong focus-visible:ring-accent h-9 w-full rounded-md border px-2 text-sm outline-none focus-visible:ring-2"
        >
          <option value="">Select a dataset...</option>
          {datasets?.map((d) => (
            <option key={d.name} value={d.name}>
              {d.name}
            </option>
          ))}
        </select>
        <div className="mt-1.5 flex flex-wrap gap-1">
          {QUICK_DATASETS.map((qd) => (
            <QuickDatasetButton
              key={qd.name}
              name={qd.name}
              label={qd.label}
              active={selectedDataset === qd.name}
              onSelect={handleQuickSelect}
            />
          ))}
        </div>
      </div>

      <div className="grid grid-cols-1 gap-2 lg:grid-cols-2">
        {FILE_SLOTS.map((slot) => {
          const required = reqs.required.includes(slot.kind);
          const optional = reqs.optional.includes(slot.kind);

          return <FileSlot key={slot.kind} config={slot} relevant={required || optional} required={required} />;
        })}
      </div>
    </div>
  );
}

function QuickDatasetButton({
  name,
  label,
  active,
  onSelect,
}: {
  name: string;
  label: string;
  active: boolean;
  onSelect: (name: string) => void;
}) {
  const handleClick = useCallback(() => onSelect(name), [name, onSelect]);

  return (
    <button
      type="button"
      onClick={handleClick}
      className={cn(
        "text-2xs rounded-sm px-1.5 py-0.5 font-mono transition-colors",
        active ? "bg-accent text-accent-fg" : "bg-surface-2 text-ink-muted hover:bg-surface-3",
      )}
    >
      {label}
    </button>
  );
}
