import { ChevronDown, ChevronRight } from "lucide-react";
import { useState, useCallback, useMemo } from "react";
import { useAppStore } from "../store/app-store";
import type { CommandName } from "../types";

interface ParamDef {
  key: string;
  label: string;
  type: "select" | "number" | "toggle" | "text";
  options?: string[];
  defaultValue: string | number | boolean;
  tooltip: string;
}

const ESSENTIAL_PARAMS: Record<CommandName, ParamDef[]> = {
  timetree: [
    {
      key: "clock_rate",
      label: "Clock rate",
      type: "text",
      defaultValue: "",
      tooltip: "Fixed clock rate (leave empty to infer)",
    },
    {
      key: "max_iter",
      label: "Max iterations",
      type: "number",
      defaultValue: 2,
      tooltip: "Maximum number of timetree iterations",
    },
    {
      key: "coalescent",
      label: "Coalescent prior",
      type: "select",
      options: ["None", "Skyline", "Constant"],
      defaultValue: "None",
      tooltip: "Coalescent tree prior",
    },
    {
      key: "reroot",
      label: "Reroot",
      type: "select",
      options: ["Least squares", "Min deviation", "Oldest", "Clock filter"],
      defaultValue: "Least squares",
      tooltip: "Rerooting strategy",
    },
  ],
  ancestral: [
    {
      key: "method_anc",
      label: "Method",
      type: "select",
      options: ["Marginal", "Parsimony"],
      defaultValue: "Marginal",
      tooltip: "Ancestral reconstruction method",
    },
    {
      key: "model_name",
      label: "Model",
      type: "select",
      options: ["Infer", "JC69", "K80", "F81", "HKY85", "T92", "TN93"],
      defaultValue: "Infer",
      tooltip: "Substitution model",
    },
    {
      key: "alphabet",
      label: "Alphabet",
      type: "select",
      options: ["Nucleotide", "Amino acid"],
      defaultValue: "Nucleotide",
      tooltip: "Sequence alphabet",
    },
  ],
  clock: [
    {
      key: "reroot",
      label: "Reroot",
      type: "select",
      options: ["Least squares", "Min deviation", "Oldest", "Clock filter", "MRCA"],
      defaultValue: "Least squares",
      tooltip: "Rerooting strategy",
    },
    {
      key: "clock_filter",
      label: "Clock filter",
      type: "number",
      defaultValue: 3.0,
      tooltip: "IQD threshold for outlier detection",
    },
  ],
  mugration: [
    {
      key: "attribute",
      label: "Attribute",
      type: "text",
      defaultValue: "country",
      tooltip: "Column name in states file to reconstruct",
    },
    { key: "iterations", label: "Iterations", type: "number", defaultValue: 5, tooltip: "GTR refinement iterations" },
  ],
  optimize: [
    {
      key: "max_iter",
      label: "Max iterations",
      type: "number",
      defaultValue: 10,
      tooltip: "Maximum optimization iterations",
    },
    { key: "dp", label: "Convergence", type: "number", defaultValue: 0.1, tooltip: "Convergence threshold" },
    {
      key: "model_name",
      label: "Model",
      type: "select",
      options: ["Infer", "JC69", "K80", "F81", "HKY85", "T92", "TN93"],
      defaultValue: "Infer",
      tooltip: "Substitution model",
    },
  ],
  prune: [
    {
      key: "prune_short",
      label: "Prune short",
      type: "text",
      defaultValue: "",
      tooltip: "Branch length threshold (e.g. 1e-6)",
    },
    {
      key: "prune_empty",
      label: "Prune empty",
      type: "toggle",
      defaultValue: false,
      tooltip: "Remove nodes with no sequence data",
    },
  ],
};

const ADVANCED_PARAMS: Record<CommandName, ParamDef[]> = {
  timetree: [
    {
      key: "branch_length_mode",
      label: "Branch lengths",
      type: "select",
      options: ["Marginal", "Input"],
      defaultValue: "Marginal",
      tooltip: "How to handle branch lengths",
    },
    {
      key: "time_marginal",
      label: "Time marginal",
      type: "select",
      options: ["Never", "Always", "Only final"],
      defaultValue: "Never",
      tooltip: "When to use marginal time estimation",
    },
    {
      key: "confidence",
      label: "Confidence intervals",
      type: "toggle",
      defaultValue: false,
      tooltip: "Compute rate-uncertainty confidence intervals",
    },
    { key: "seed", label: "Random seed", type: "text", defaultValue: "", tooltip: "Seed for reproducibility" },
  ],
  ancestral: [
    { key: "dense", label: "Dense mode", type: "toggle", defaultValue: false, tooltip: "Force dense representation" },
    {
      key: "gap_fill",
      label: "Gap fill",
      type: "select",
      options: ["Only terminal", "All", "None"],
      defaultValue: "Only terminal",
      tooltip: "Gap filling strategy",
    },
    {
      key: "reconstruct_tip_states",
      label: "Reconstruct tips",
      type: "toggle",
      defaultValue: false,
      tooltip: "Reconstruct tip states",
    },
    { key: "seed", label: "Random seed", type: "text", defaultValue: "", tooltip: "Seed for reproducibility" },
  ],
  clock: [
    {
      key: "branch_length_mode",
      label: "Branch lengths",
      type: "select",
      options: ["Marginal", "Input"],
      defaultValue: "Marginal",
      tooltip: "How to handle branch lengths",
    },
    {
      key: "keep_root",
      label: "Keep root",
      type: "toggle",
      defaultValue: false,
      tooltip: "Keep the original root position",
    },
    {
      key: "covariation",
      label: "Covariation",
      type: "toggle",
      defaultValue: false,
      tooltip: "Account for covariation in rate estimation",
    },
    { key: "seed", label: "Random seed", type: "text", defaultValue: "", tooltip: "Seed for reproducibility" },
  ],
  mugration: [
    {
      key: "pc",
      label: "Pseudo counts",
      type: "text",
      defaultValue: "",
      tooltip: "Pseudo counts (leave empty to compute)",
    },
    {
      key: "missing_data",
      label: "Missing symbol",
      type: "text",
      defaultValue: "?",
      tooltip: "Symbol for missing data",
    },
    {
      key: "sampling_bias_correction",
      label: "Sampling bias",
      type: "text",
      defaultValue: "",
      tooltip: "Sampling bias correction factor",
    },
  ],
  optimize: [
    {
      key: "damping",
      label: "Damping",
      type: "number",
      defaultValue: 0.75,
      tooltip: "Damping factor for branch optimization",
    },
    {
      key: "opt_method",
      label: "Optimization",
      type: "select",
      options: ["Brent sqrt", "Brent", "Brent log", "Newton", "Newton sqrt", "Newton log"],
      defaultValue: "Brent sqrt",
      tooltip: "Branch length optimization method",
    },
    {
      key: "no_indels",
      label: "Ignore indels",
      type: "toggle",
      defaultValue: false,
      tooltip: "Ignore insertions and deletions",
    },
  ],
  prune: [
    {
      key: "merge_shared_mutations",
      label: "Merge shared",
      type: "toggle",
      defaultValue: false,
      tooltip: "Merge shared mutations",
    },
  ],
};

export function ParamForm() {
  const activeCommand = useAppStore((s) => s.activeCommand);
  const [showAdvanced, setShowAdvanced] = useState(false);

  const essentials = useMemo(() => ESSENTIAL_PARAMS[activeCommand], [activeCommand]);
  const advanced = useMemo(() => ADVANCED_PARAMS[activeCommand], [activeCommand]);

  const toggleAdvanced = useCallback(() => setShowAdvanced((s) => !s), []);

  return (
    <div className="space-y-3">
      <h3 className="text-xs font-semibold uppercase tracking-wider text-gray-500 dark:text-gray-400">Parameters</h3>

      <div className="space-y-2">
        {essentials.map((param) => (
          <ParamField key={param.key} param={param} />
        ))}
      </div>

      {advanced.length > 0 && (
        <div>
          <button
            type="button"
            onClick={toggleAdvanced}
            className="flex items-center gap-1 text-xs font-medium text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-300"
          >
            {showAdvanced ? <ChevronDown size={14} /> : <ChevronRight size={14} />}
            Advanced
          </button>
          {showAdvanced && (
            <div className="mt-2 space-y-2 border-l-2 border-gray-200 pl-3 dark:border-gray-700">
              {advanced.map((param) => (
                <ParamField key={param.key} param={param} />
              ))}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

function ParamField({ param }: { param: ParamDef }) {
  const [value, setValue] = useState(param.defaultValue);

  return (
    <div className="flex items-center gap-3" title={param.tooltip}>
      <label className="w-28 shrink-0 text-xs text-gray-600 dark:text-gray-400">{param.label}</label>

      {param.type === "select" && (
        <select
          value={String(value)}
          onChange={(e) => setValue(e.target.value)}
          className="flex-1 rounded-md border border-gray-200 bg-white px-2 py-1 text-xs text-gray-700 dark:border-gray-600 dark:bg-gray-800 dark:text-gray-300"
        >
          {param.options?.map((opt) => (
            <option key={opt} value={opt}>
              {opt}
            </option>
          ))}
        </select>
      )}

      {param.type === "number" && (
        <input
          type="number"
          value={String(value)}
          onChange={(e) => setValue(Number(e.target.value))}
          step="any"
          className="flex-1 rounded-md border border-gray-200 bg-white px-2 py-1 text-xs text-gray-700 dark:border-gray-600 dark:bg-gray-800 dark:text-gray-300"
        />
      )}

      {param.type === "text" && (
        <input
          type="text"
          value={String(value)}
          onChange={(e) => setValue(e.target.value)}
          placeholder={param.tooltip}
          className="flex-1 rounded-md border border-gray-200 bg-white px-2 py-1 text-xs text-gray-700 placeholder:text-gray-400 dark:border-gray-600 dark:bg-gray-800 dark:text-gray-300 dark:placeholder:text-gray-600"
        />
      )}

      {param.type === "toggle" && (
        <button
          type="button"
          onClick={() => setValue(!value)}
          className={`relative h-5 w-9 rounded-full transition-colors ${
            value ? "bg-[var(--color-accent)]" : "bg-gray-300 dark:bg-gray-600"
          }`}
        >
          <span
            className={`absolute top-0.5 left-0.5 h-4 w-4 rounded-full bg-white shadow transition-transform ${
              value ? "translate-x-4" : ""
            }`}
          />
        </button>
      )}
    </div>
  );
}
