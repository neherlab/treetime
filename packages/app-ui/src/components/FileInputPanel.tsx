import { FileSlot } from "./FileSlot";
import { useAppStore } from "../store/app-store";
import type { CommandName, FileSlotKind } from "../types";
import { FILE_SLOTS } from "../types";

const COMMAND_FILE_REQUIREMENTS: Record<CommandName, { required: FileSlotKind[]; optional: FileSlotKind[] }> = {
  timetree: {
    required: ["tree", "alignment", "dates"],
    optional: ["vcfReference"],
  },
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
};

export function FileInputPanel() {
  const activeCommand = useAppStore((s) => s.activeCommand);
  const reqs = COMMAND_FILE_REQUIREMENTS[activeCommand];

  return (
    <div className="space-y-2">
      <h3 className="text-xs font-semibold uppercase tracking-wider text-gray-500 dark:text-gray-400">Input files</h3>
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
