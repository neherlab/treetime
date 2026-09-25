import type { AppCommand, InputKind } from "@neherlab/app-contracts";

interface CommandInfo {
  command: AppCommand;
  label: string;
  verb: string;
  description: string;
}

export const APP_COMMANDS: readonly AppCommand[] = ["timetree", "clock", "ancestral", "mugration", "optimize", "prune"];

export const COMMAND_INFO: Record<AppCommand, CommandInfo> = {
  timetree: {
    command: "timetree",
    label: "Time tree",
    verb: "Run time tree",
    description: "Date the ancestors and estimate the clock rate from sampling dates.",
  },
  clock: {
    command: "clock",
    label: "Clock signal",
    verb: "Run clock check",
    description: "Root-to-tip regression: clock rate, temporal signal and outliers.",
  },
  ancestral: {
    command: "ancestral",
    label: "Ancestral sequences",
    verb: "Reconstruct sequences",
    description: "Infer ancestral sequences and the mutations on each branch.",
  },
  mugration: {
    command: "mugration",
    label: "Discrete traits",
    verb: "Reconstruct traits",
    description: "Ancestral states of a metadata column, such as country or host.",
  },
  optimize: {
    command: "optimize",
    label: "Branch lengths",
    verb: "Optimize branch lengths",
    description: "Maximum-likelihood branch lengths on a fixed topology.",
  },
  prune: {
    command: "prune",
    label: "Prune tree",
    verb: "Prune tree",
    description: "Collapse short or empty branches and remove listed samples.",
  },
};

export const INPUT_SLOT_INFO: Record<InputKind, { label: string; hint: string; extensions: string[] }> = {
  tree: { label: "Tree", hint: "Newick or Nexus", extensions: ["nwk", "newick", "nex", "nexus", "tree", "tre"] },
  alignment: {
    label: "Alignment",
    hint: "FASTA, aligned",
    extensions: ["fasta", "fa", "fas", "aln", "xz", "gz", "zst", "bz2"],
  },
  metadata: { label: "Metadata", hint: "TSV or CSV with names and dates", extensions: ["tsv", "csv", "txt"] },
};

export const MAIN_SETTING_KEYS: Record<AppCommand, readonly string[]> = {
  timetree: [
    "clock_rate",
    "clock_std_dev",
    "confidence",
    "covariation",
    "coalescent",
    "coalescent_opt",
    "coalescent_skyline",
    "skyline_n_points",
    "skyline_stiffness",
    "reroot",
    "keep_root",
    "clock_filter",
    "relax",
    "keep_polytomies",
    "model",
    "model_params",
    "max_iter",
  ],
  clock: ["reroot", "keep_root", "clock_filter", "covariation", "allow_negative_rate", "metadata_id_columns"],
  ancestral: ["method_anc", "model", "model_params", "gap_fill", "reconstruct_tip_states"],
  mugration: ["attribute", "weights", "pc", "missing_data", "sampling_bias_correction"],
  optimize: ["opt_method", "reroot", "divergence_units", "no_indels", "max_iter"],
  prune: ["prune_short", "prune_empty", "merge_shared_mutations", "prune_nodes_list"],
};

export function isAppCommand(value: string | undefined): value is AppCommand {
  return APP_COMMANDS.some((command) => command === value);
}
