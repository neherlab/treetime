import { RotateCcw } from "lucide-react";
import { useCallback } from "react";
import { FileSlot } from "./FileSlot";
import { useDatasets } from "../hooks";
import { useAppStore } from "../store/app-store";
import type { CommandName, FileSlotKind } from "../types";
import { FILE_SLOTS } from "../types";

const COMMAND_FILE_REQUIREMENTS: Record<CommandName, { required: FileSlotKind[]; optional: FileSlotKind[] }> = {
  ancestral: {
    required: ["tree", "alignment"],
    optional: ["vcfReference"],
  },
  clock: {
    required: ["tree", "dates"],
    optional: ["alignment", "vcfReference"],
  },
  mugration: {
    required: ["tree", "states"],
    optional: ["weights"],
  },
  optimize: {
    required: ["tree", "alignment"],
    optional: [],
  },
  prune: {
    required: ["tree"],
    optional: ["alignment"],
  },
  timetree: {
    required: ["tree", "alignment", "dates"],
    optional: ["vcfReference"],
  },
};

export function FileInputPanel() {
  const activeCommand = useAppStore((s) => s.activeCommand);
  const selectedDataset = useAppStore((s) => s.selectedDataset);
  const setSelectedDataset = useAppStore((s) => s.setSelectedDataset);
  const setFile = useAppStore((s) => s.setFile);
  const resetForm = useAppStore((s) => s.resetForm);
  const reqs = COMMAND_FILE_REQUIREMENTS[activeCommand];
  const { data: datasets } = useDatasets();

  const handleDatasetChange = useCallback(
    (e: React.ChangeEvent<HTMLSelectElement>) => {
      const dataset = e.target.value;
      setSelectedDataset(dataset);
      if (dataset) {
        const info = datasets?.find((d) => d.name === dataset);
        const files = info?.files ?? [];
        if (files.includes("tree.nwk")) {
          setFile("tree", { name: `${dataset}/tree.nwk`, size: 0 });
        }
        if (files.includes("aln.fasta.xz")) {
          setFile("alignment", { name: `${dataset}/aln.fasta.xz`, size: 0 });
        }
        if (files.includes("metadata.tsv")) {
          setFile("dates", { name: `${dataset}/metadata.tsv`, size: 0 });
          setFile("states", { name: `${dataset}/metadata.tsv`, size: 0 });
        }
      }
    },
    [setSelectedDataset, setFile, datasets],
  );

  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <h3 className="text-xs font-semibold uppercase tracking-wider text-gray-500 dark:text-gray-400">Input files</h3>
        <button
          type="button"
          onClick={resetForm}
          className="flex items-center gap-1 rounded-md px-2 py-0.5 text-xs text-gray-500 hover:bg-gray-100 hover:text-gray-700 dark:text-gray-400 dark:hover:bg-gray-800 dark:hover:text-gray-200"
        >
          <RotateCcw size={12} />
          Reset
        </button>
      </div>

      <div>
        <label className="mb-1 block text-xs text-gray-500 dark:text-gray-400">Example dataset</label>
        <select
          value={selectedDataset}
          onChange={handleDatasetChange}
          className="w-full rounded-md border border-gray-200 bg-white px-2 py-1.5 text-xs text-gray-700 dark:border-gray-600 dark:bg-gray-800 dark:text-gray-300"
        >
          <option value="">Select a dataset...</option>
          {datasets?.map((d) => (
            <option key={d.name} value={d.name}>
              {d.name}
            </option>
          ))}
        </select>
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
